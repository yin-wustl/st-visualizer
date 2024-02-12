import { colTypes, importStateType } from "../imports/api/constants";
import * as fs from "fs/promises";
export const ipcConstants = {
  doCalculation: "doCalculation",
  getFile: "getFile",
} as const;

export const ipcRenderConstants = {
  saveFile: "saveFile",
  openFile: "openFile",
};

import path from "path";
import { execFile } from "child_process";
import _ from "lodash";
import { app } from "electron";

export const ipcHandlers = {
  [ipcConstants.getFile]: async (
    _event: Electron.IpcMainInvokeEvent,
    { path }: { path: string }
  ) => {
    const string = await fs.readFile(path, { encoding: "utf-8" });
    return string;
  },
  [ipcConstants.doCalculation]: async (
    _event: Electron.IpcMainInvokeEvent,
    {
      transforms,
      importState,
    }: {
      transforms: { x: number; y: number }[][];
      importState: importStateType;
    }
  ) => {
    const staticsPath = app.isPackaged
      ? path.resolve(process.resourcesPath)
      : path.resolve(__dirname, "../../imports/static/");

    const exePath = process.platform === "win32"
      ? path.resolve(staticsPath, "./st-visualizer.exe")
      : path.resolve(staticsPath, "./st-visualizer");

    const split = transforms.map(
      (slice) =>
        [slice.map((v) => v.x), slice.map((v) => v.y)] as [number[], number[]]
    );

    const baseLayer = split.at(0);
    if (!baseLayer) return;
    const rest = split?.slice(1);
    const alignments =
      " , , , , , ,\n" +
      rest
        ?.flatMap((r) => [r, baseLayer])
        .flatMap((rowPair) => rowPair)
        .map((row) => " , , ," + row.join(","))
        .join("\n");

    const alignmentPath = path.resolve(__dirname, "./alignments.csv");

    await fs.writeFile(alignmentPath, alignments, {
      flag: "w",
    });
    const outputPath = path.resolve(__dirname, "./output.json");
    const sliceNames = _.compact(
      _.uniq(
        importState.tsvData.map((row, i) =>
          i > 0 ? row[importState[colTypes.slice]] : undefined
        )
      )
    );

    try {
      const configOptions = {
        fileName: importState.fileName,
        alignmentFile: alignmentPath,
        target: outputPath,
        shrink: importState.shrink,
        sliceNames: importState.sliceOrder.map((i) => sliceNames[i]),
        featureCols: importState[colTypes.feature],
        sliceIndex: importState[colTypes.slice],
        tissueIndex: importState[colTypes.tissue],
        rowIndex: importState[colTypes.row],
        colIndex: importState[colTypes.column],
        clusterIndex: importState[colTypes.cluster],
        zDistance: importState.z_distance,
        resultExport: true,
        objExport: false,
        featureObj: null,
        clusterObj: null,
      };
      console.log(JSON.stringify(configOptions));
      const runResults = await new Promise<string>((res, rej) => {
        execFile(
          exePath,
          ["1", JSON.stringify(configOptions)],
          {},
          (err, data) => {
            if (err) rej(err);
            res(data);
          }
        );
      });
      console.log(runResults);
      // console.log("___COMPLETE___");
      // console.log("Data output to ", outputPath);

      return outputPath;
    } catch (e) {
      console.error("Running error. Please try again: ", e);
    }
  },
} as const;
