import json

def parseParameters(jsonFile):
    with open(jsonFile) as f:
        return json.load(f)