import argparse

import json



import src.download_L4 as dw



parser = argparse.ArgumentParser()

parser.add_argument("--json", help="insert parameters oilspill json format")    

args = parser.parse_args()



with open(args.json) as f:

    jdata = json.load(f)





dw.download_cmems_data(jdata)
