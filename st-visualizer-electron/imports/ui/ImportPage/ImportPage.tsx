import {
  Button,
  FormControl,
  Grid,
  InputLabel,
  MenuItem,
  Select,
  Stack,
  TextField,
  Typography,
} from "@mui/material";
import React, { useState } from "react";
import FileOpenIcon from '@mui/icons-material/FileOpen';
import { getValue, importStateType, colTypes } from "../../api/constants";
import _ from "lodash";

export const ImportPage = ({
  importState,
  setImportState,
}: {
  importState: importStateType;
  setImportState: React.Dispatch<React.SetStateAction<importStateType>>;
}): JSX.Element => {
  const tsvLabels = importState.tsvData.at(0) || [];

  const nCols = tsvLabels.length;

  const changeSelection = (i: number, selection: number) => {
    const currentValue = getValue(importState, i);
    if (currentValue === colTypes.feature) {
      setImportState((s) => ({
        ...s,
        [colTypes.feature]: s[colTypes.feature].filter((v) => v !== i),
      }));
    } else {
      setImportState((s) => ({ ...s, [currentValue]: -1 }));
    }
    if (selection === colTypes.feature) {
      setImportState((s) => ({
        ...s,
        [colTypes.feature]: [...s[colTypes.feature], i],
      }));
    } else {
      setImportState((s) => ({ ...s, [selection]: i }));
    }
  };

  const columns = (
    <Grid item container>
      <Typography variant="h4" style={{ paddingTop: "15px", paddingBottom: "15px" }}>Column Management</Typography>
      <Grid item container spacing={1}>
        {[...new Array(nCols)].map((e, i) => (
          <Grid item key={i}>
            <FormControl>
              <InputLabel>{tsvLabels[i]}</InputLabel>
              <Select
                label={tsvLabels[i]}
                style={{ minWidth: 130 }}
                value={getValue(importState, i)}
                onChange={(e) => {
                  changeSelection(i, e.target.value as colTypes);
                  if (e.target.value === colTypes.slice) {
                    const slicesRow = importState.tsvData.map(
                      (row) => row[i]
                    );
                    // console.log(slicesRow)
                    const numSlices = _.compact(_.uniq(slicesRow.slice(1))).length;
                    setImportState((s) => ({
                      ...s,
                      numSlices: numSlices,
                    }));
                  }
                }}
              >
                <MenuItem value={colTypes.none}>
                  <em>None</em>
                </MenuItem>
                <MenuItem value={colTypes.slice}>Slice</MenuItem>
                <MenuItem value={colTypes.tissue}>Tissue</MenuItem>
                <MenuItem value={colTypes.row}>Row</MenuItem>
                <MenuItem value={colTypes.column}>Column</MenuItem>
                <MenuItem value={colTypes.cluster}>Cluster</MenuItem>
                <MenuItem value={colTypes.feature}>Feature</MenuItem>
              </Select>
            </FormControl>
          </Grid>
        ))}
      </Grid>
    </Grid>
  );
  const metadata = (
    <Grid item container>
      <Typography variant="h4" style={{ paddingTop: "15px", paddingBottom: "15px" }}>Configurations</Typography>
      <Grid item container spacing={1}>
        <Grid item>
          <TextField
            id="outlined-number"
            type="number"
            label={"Shrink (microns)"}
            value={importState.shrink}
            onChange={(e) =>
              setImportState((s) => ({ ...s, shrink: parseFloat(e.target.value) }))
            }
            InputLabelProps={{
              shrink: true,
            }}
          />
        </Grid>
        <Grid item>
          <TextField
            //Set of widths rather than a single constant one
            id="outlined-number"
            type="number"
            label={"Slice Width (Microns)"}
            value={importState.z_distance}
            onChange={(e) =>
              setImportState((s) => ({
                ...s,
                z_distance: parseInt(e.target.value),
              }))
            }
            InputLabelProps={{
              shrink: true,
            }}
          />
        </Grid>
      </Grid>
    </Grid>
  );

  return (
    <Stack style={{ padding: 20 }}>
      <Button variant="contained" component="label" startIcon={<FileOpenIcon />} style={{textTransform: 'none'}}>
        Data File {importState.fileName ? importState.fileName : null}
        <input
          hidden
          multiple={false}
          accept=".tsv"
          type="file"
          onChange={async (e) => {
            const file = e.target?.files?.[0];
            if (!file) return;
            const tsvPath = file.path;
            const rawText = await file.text();
            const tsvData = rawText.split("\n").filter((l) => l.length > 1).map((l) => l.split("\t"));
            const slicesRow = tsvData.map(
              (row) => row[importState[colTypes.slice]]
            );
            const numSlices = _.compact(_.uniq(slicesRow.slice(1))).length;
            setImportState((s) => ({
              ...s,
              numSlices: numSlices,
              tsvData: tsvData,
              fileName: tsvPath,
              [colTypes.feature]: [...new Array(tsvData[0].length - 6)].map(
                (_, i) => i + 6
              ),
            }));
          }}
        />
      </Button>
      {importState.fileName ? columns : null}
      {importState.fileName ? metadata : null}
    </Stack>
  );
};
