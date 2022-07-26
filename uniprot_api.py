import re
import time
import json
import zlib
from xml.etree import ElementTree
from urllib.parse import urlparse, parse_qs, urlencode
import requests
from requests.adapters import HTTPAdapter, Retry
import csv
import pandas as pd

POLLING_INTERVAL = 3
API_URL = "https://rest.uniprot.org"


retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))


def check_response(response):
    try:
        response.raise_for_status()
    except requests.HTTPError:
        print(response.json())
        raise


def submit_id_mapping(from_db, to_db, ids):
    request = requests.post(
        f"{API_URL}/idmapping/run",
        data={"from": from_db, "to": to_db, "ids": ",".join(ids)},
    )
    check_response(request)
    return request.json()["jobId"]


def get_next_link(headers):
    re_next_link = re.compile(r'<(.+)>; rel="next"')
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)


def check_id_mapping_results_ready(job_id):
    while True:
        request = session.get(f"{API_URL}/idmapping/status/{job_id}")
        check_response(request)
        j = request.json()
        if "jobStatus" in j:
            if j["jobStatus"] == "RUNNING":
                print(f"Retrying in {POLLING_INTERVAL}s")
                time.sleep(POLLING_INTERVAL)
            else:
                raise Exception(j["jobStatus"])
        else:
            return bool(j["results"] or j["failedIds"])


def get_batch(batch_response, file_format, compressed):
    batch_url = get_next_link(batch_response.headers)
    while batch_url:
        batch_response = session.get(batch_url)
        batch_response.raise_for_status()
        yield decode_results(batch_response, file_format, compressed)
        batch_url = get_next_link(batch_response.headers)


def combine_batches(all_results, batch_results, file_format):
    if file_format == "json":
        for key in ("results", "failedIds"):
            if key in batch_results and batch_results[key]:
                all_results[key] += batch_results[key]
    elif file_format == "tsv":
        return all_results + batch_results[1:]
    else:
        return all_results + batch_results
    return all_results


def get_id_mapping_results_link(job_id):
    url = f"{API_URL}/idmapping/details/{job_id}"
    request = session.get(url)
    check_response(request)
    return request.json()["redirectURL"]


def decode_results(response, file_format, compressed):
    if compressed:
        decompressed = zlib.decompress(response.content, 16 + zlib.MAX_WBITS)
        if file_format == "json":
            j = json.loads(decompressed.decode("utf-8"))
            return j
        elif file_format == "tsv":
            return [line for line in decompressed.decode("utf-8").split("\n") if line]
        elif file_format == "xlsx":
            return [decompressed]
        elif file_format == "xml":
            return [decompressed.decode("utf-8")]
        else:
            return decompressed.decode("utf-8")
    elif file_format == "json":
        return response.json()
    elif file_format == "tsv":
        return [line for line in response.text.split("\n") if line]
    elif file_format == "xlsx":
        return [response.content]
    elif file_format == "xml":
        return [response.text]
    return response.text


def get_xml_namespace(element):
    m = re.match(r"\{(.*)\}", element.tag)
    return m.groups()[0] if m else ""


def merge_xml_results(xml_results):
    merged_root = ElementTree.fromstring(xml_results[0])
    for result in xml_results[1:]:
        root = ElementTree.fromstring(result)
        for child in root.findall("{http://uniprot.org/uniprot}entry"):
            merged_root.insert(-1, child)
    ElementTree.register_namespace("", get_xml_namespace(merged_root[0]))
    return ElementTree.tostring(merged_root, encoding="utf-8", xml_declaration=True)


def print_progress_batches(batch_index, size, total):
    n_fetched = min((batch_index + 1) * size, total)
    print(f"Fetched: {n_fetched} / {total}")


def get_id_mapping_results_search(url):
    parsed = urlparse(url)
    query = parse_qs(parsed.query)
    file_format = query["format"][0] if "format" in query else "json"
    if "size" in query:
        size = int(query["size"][0])
    else:
        size = 500
        query["size"] = size
    compressed = (
        query["compressed"][0].lower() == "true" if "compressed" in query else False
    )
    parsed = parsed._replace(query=urlencode(query, doseq=True))
    url = parsed.geturl()
    request = session.get(url)
    check_response(request)
    results = decode_results(request, file_format, compressed)
    total = int(request.headers["x-total-results"])
    print_progress_batches(0, size, total)
    for i, batch in enumerate(get_batch(request, file_format, compressed), 1):
        results = combine_batches(results, batch, file_format)
        print_progress_batches(i, size, total)
    if file_format == "xml":
        return merge_xml_results(results)
    return results


def get_id_mapping_results_stream(url):
    if "/stream/" not in url:
        url = url.replace("/results/", "/results/stream/")
    request = session.get(url)
    check_response(request)
    parsed = urlparse(url)
    query = parse_qs(parsed.query)
    file_format = query["format"][0] if "format" in query else "json"
    compressed = (
        query["compressed"][0].lower() == "true" if "compressed" in query else False
    )
    return decode_results(request, file_format, compressed)


def get_data_frame_from_tsv_results(tsv_results):
    reader = csv.DictReader(tsv_results, delimiter="\t", quotechar='"')
    return pd.DataFrame(list(reader))


def api_map_STRING_GeneName(ppi):
    string_ids = ppi.iloc[:, 0].drop_duplicates().tolist()
    job_id = submit_id_mapping(from_db="STRING", to_db="UniProtKB", ids=string_ids)
    if check_id_mapping_results_ready(job_id):
        link = get_id_mapping_results_link(job_id)
        results = get_id_mapping_results_search(link)
        mapping_dict = dict() # mapping_dict[<string_id>] = <UniProtKB_geneName>
        for entry in results['results']:
            try:
                    # mapping_dict[entry['from']] = entry['to']['genes'][0]['geneName']['value']
                if entry['to']['entryType'] == 'UniProtKB unreviewed (TrEMBL)':
                    mapping_dict[entry['from']] = entry['to']['primaryAccession'] # UniProtKB unreviewed ID, not in results['failedIds']
                    print(mapping_dict[entry['from']], entry['to']['entryType'])
            except KeyError:
                try:
                    mapping_dict[entry['from']] = entry['to']['genes'][0]['orderedLocusNames'][0]['value']
                except KeyError:
                    continue
    fh = open('/home/ahphan/RotationData/Friedberg/bar/unknown-strings.txt', "w")
    fh2 = open('/home/ahphan/RotationData/Friedberg/bar/trembl-strings-names.txt', "w")
    for string_id in mapping_dict.keys():
        fh.write(string_id + "\n")
        fh2.write(string_id+": "+ mapping_dict[string_id]+"\n")
    fh.close
    fh2.close
    return mapping_dict, results


def api_map_STRING_UniProt(ppi):
    string_ids = ppi.iloc[:, 0].drop_duplicates().tolist() # STRING IDs messed up because UniProt API have missing fields
    # On string-db.org, it points to different protein than when mapped with API
    job_id = submit_id_mapping(from_db="STRING", to_db="UniProtKB", ids=string_ids)
    if check_id_mapping_results_ready(job_id):
        link = get_id_mapping_results_link(job_id)
        results = get_id_mapping_results_search(link)
        mapping_dict = dict() # mapping_dict[<string_id>] = <UniProtKB_ID>
        mapping_dict_unre = dict()
        for entry in results['results']:
            try:
                    # mapping_dict[entry['from']] = entry['to']['genes'][0]['geneName']['value']
                if entry['to']['entryType'] == 'UniProtKB unreviewed (TrEMBL)':
                    mapping_dict_unre[entry['from']] = entry['to']['primaryAccession'] # UniProtKB unreviewed ID, not in results['failedIds']
                    print(mapping_dict_unre[entry['from']], entry['to']['entryType'])
                else:
                    mapping_dict[entry['from']] = entry['to']['primaryAccession']
            except KeyError:
                try:
                    mapping_dict[entry['from']] = entry['to']['genes'][0]['orderedLocusNames'][0]['value']
                except KeyError:
                    continue
    fh = open('unknown-strings.txt', "w")
    fh2 = open('trembl-strings-names.txt', "w")
    for string_id in mapping_dict_unre.keys():
        fh.write(string_id + "\n")
        fh2.write(string_id+":"+ mapping_dict_unre[string_id]+"\n")
    fh.close
    fh2.close
    return mapping_dict, mapping_dict_unre

# ppi = pd.read_csv('/home/ahphan/Downloads/9606.protein.links.coexp400.v11.5.txt', header=0,
#                       delim_whitespace=True)
# mapping_dict = api_map_STRING_GeneName(ppi)

# mapping_table = pd.DataFrame.from_dict(mapping_dict, orient='index')
# print(mapping_table.shape)
# mapping_na = mapping_table[mapping_table.isna().any(axis=1)]
# print(mapping_na.head())
# print(mapping_na.shape)

    # Equivalently using the stream endpoint which is more demanding
    # on the API and so is less stable:
    # results = get_id_mapping_results_stream(link)

# print(results)

# print(mapping_dict['511145.b0360'])
# results2 = get_id_mapping_results_stream("https://rest.uniprot.org/idmapping/uniprotkb/results/stream/604970e6779a6288c1d7dd9a932bedc7eee056b8?compressed=true&fields=accession%2Cgene_names&format=tsv")
# print(results2)