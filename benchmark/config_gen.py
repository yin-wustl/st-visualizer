import json

def config_gen(size: str, path: str, alignment: str):
    output = {
        "fileName": f"{path}data_{size}.tsv",
        "alignmentFile": alignment,
        "target": f"{path}output_{size}.json",
        "shrink": 0,
        "sliceNames": [
            "slice_0",
            "slice_1",
            "slice_2",
            "slice_3",
            "slice_4",
            "slice_5",
            "slice_6",
            "slice_7",
            "slice_8",
            "slice_9",
        ],
        "featureCols": [5, 6, 7, 8, 9, 10, 11, 12, 13, 14],
        "sliceIndex": 0,
        "tissueIndex": 1,
        "rowIndex": 16,
        "colIndex": 15,
        "clusterIndex": 2,
        "zDistance": 10,
        "resultExport": False,
        "objExport": False,
        "featureObj": f"{path}feature_{size}/",
        "clusterObj": f"{path}cluster_{size}/",
    }

    # write to file
    with open(f"{path}config_{size}.json", "w") as file:
        json.dump(output, file)
