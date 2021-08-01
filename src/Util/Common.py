import json
import jsonlines


def jsonlToJsonConverter(readFilePath, writeFilePath):
    """ file.jsonl -> file.json"""
    with jsonlines.open(readFilePath, "r") as rfd:
        with open(writeFilePath, "w", encoding="utf-8") as wfd:
            for data in rfd:
                json.dump(data, wfd, indent=4, ensure_ascii=False)
