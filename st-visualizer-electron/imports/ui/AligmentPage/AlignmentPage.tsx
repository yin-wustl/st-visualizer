import React, { useEffect, useRef, useState } from "react";
import {
  Stack,
  Button,
  Grid,
  IconButton,
  Tooltip,
  Slider,
  Input,
  Typography,
  Box,
  Paper,
} from "@mui/material";
import { ArrowLeft, ArrowRight } from "@mui/icons-material";
import ImageIcon from "@mui/icons-material/Image";
import OpacityIcon from "@mui/icons-material/Opacity";
import SwapHorizIcon from "@mui/icons-material/SwapHoriz";
import SwapVertIcon from "@mui/icons-material/SwapVert";
import LoopIcon from "@mui/icons-material/Loop";
import ClearIcon from "@mui/icons-material/Clear";
import ReplayIcon from "@mui/icons-material/Replay";
import { importStateType, colTypes } from "../../api/constants";
import _ from "lodash";

export type transformType = {
  file: File;
  alignments: {
    x: number;
    y: number;
    rotZ: number;
  };
};

const getImg = (src: string) =>
  new Promise<HTMLImageElement>((res) => {
    const img = new Image();
    img.src = src;
    img.onload = () => {
      res(img);
    };
  });

export const calcTransforms = (currentImages: transformType[]) => {
  const startRadius = 100;
  const startPointLocations = [
    { x: 0, y: 0 },
    { x: startRadius, y: startRadius },
    { x: -startRadius, y: -startRadius },
    { x: -startRadius, y: startRadius },
    { x: startRadius, y: -startRadius },
  ];
  const rotate = ({ x, y, angle }: { x: number; y: number; angle: number }) => {
    const radAngle = (angle * Math.PI) / 180;
    return {
      x: x * Math.cos(radAngle) + y * Math.sin(radAngle),
      y: y * Math.cos(radAngle) - x * Math.sin(radAngle),
    };
  };
  const transforms = currentImages.reduceRight((arr, { alignments }) => {
    const appendedArray = [startPointLocations, ...arr];
    //Rotate then transform
    const transformed = appendedArray.map((prev) =>
      prev
        .map((rotated) => ({
          x: rotated.x - alignments.x,
          y: rotated.y - alignments.y,
        }))
        .map((p) =>
          rotate({
            x: p.x,
            y: p.y,
            angle: alignments.rotZ,
          })
        )
    );
    return transformed;
  }, [] as { x: number; y: number }[][]);
  return transforms.map((slice) => slice.map((pt) => ({ x: pt.x, y: pt.y })));
};

export const AlignmentPage = ({
  importState,
  setImportState,
  currentImages,
  setCurrentImages,
}: {
  importState: importStateType;
  setImportState: React.Dispatch<React.SetStateAction<importStateType>>;
  currentImages: transformType[];
  setCurrentImages: React.Dispatch<React.SetStateAction<transformType[]>>;
}): JSX.Element => {
  const sliceNames = _.compact(
    _.uniq(
      importState.tsvData.map((row, i) =>
        i > 0 ? row[importState[colTypes.slice]] : undefined
      )
    )
  );
  const [selectedImg, setSelectedImg] = useState<number | null>(null);
  const [imgUrls, setImgUrls] = useState<string[]>([]);
  const transforms = calcTransforms(currentImages);

  //Slice to remove the title and the compact to remove undefined
  const numSlices = importState.numSlices;

  // console.log("CURRENT Transforms", transforms, sliceNames);
  const xCol = importState[colTypes.column];
  const yCol = importState[colTypes.row];
  const currentVals = importState.tsvData
    .map((row, i) =>
      i === 0
        ? undefined
        : {
            x: parseFloat(row[xCol]),
            y: parseFloat(row[yCol]),
            slice: importState.tsvData[importState[colTypes.slice]],
          }
    )
    .filter((v) => !!v);
  // console.log(currentVals);
  useEffect(() => {
    const func = async () => {
      // console.log("RUNNING");
      const urls = await Promise.all(
        currentImages.map(
          (file) =>
            new Promise<string>((res) => {
              const reader = new FileReader();
              reader.readAsDataURL(file.file);
              reader.onloadend = function () {
                res(reader.result as string);
              };
            })
        )
      );
      setImgUrls(urls);
    };
    func();
  }, [currentImages]);

  const imageOrdering = (
    <Stack style={{ paddingTop: "15px", paddingBottom: "15px" }}>
      <Stack direction="row">
        {imgUrls?.map((img, i) => {
          return (
            <Stack
              direction="row"
              key={img}
              style={{
                borderBottom:
                  selectedImg !== null &&
                  (i === selectedImg || i === selectedImg - 1)
                    ? "solid 5px black"
                    : undefined,
              }}
            >
              <Grid item>
                <IconButton
                  onClick={() => {
                    if (i === 0) return;
                    const temp = [...currentImages];
                    [temp[i], temp[i - 1]] = [temp[i - 1], temp[i]];
                    setCurrentImages(temp);
                  }}
                >
                  <ArrowLeft />
                </IconButton>
              </Grid>
              <Tooltip key={img} title={currentImages[i]?.file.name}>
                <Grid
                  item
                  sx={{ width: 100 }}
                  onClick={() => i && setSelectedImg(i)}
                >
                  <img
                    alt={currentImages[i]?.file.name}
                    src={img}
                    style={{ maxWidth: "100%" }}
                  />
                </Grid>
              </Tooltip>
              <Grid item>
                <IconButton
                  onClick={() => {
                    if (i === currentImages.length - 1) return;
                    const temp = [...currentImages];
                    [temp[i], temp[i + 1]] = [temp[i + 1], temp[i]];
                    setCurrentImages(temp);
                  }}
                >
                  <ArrowRight />
                </IconButton>
              </Grid>
            </Stack>
          );
        })}
      </Stack>
      <Stack direction="row">
        {importState.sliceOrder?.map((nameIndex, i) => {
          return (
            <Stack direction="row" key={nameIndex}>
              <Grid item>
                <IconButton
                  onClick={() => {
                    if (i === 0) return;
                    const temp = [...importState.sliceOrder];
                    [temp[i], temp[i - 1]] = [temp[i - 1], temp[i]];
                    setImportState((s) => ({ ...s, sliceOrder: temp }));
                  }}
                >
                  <ArrowLeft />
                </IconButton>
              </Grid>
              <Grid item sx={{ width: 100 }}>
                {sliceNames[nameIndex]}
              </Grid>
              <Grid item>
                <IconButton
                  onClick={() => {
                    if (i === currentImages.length - 1) return;
                    const temp = [...importState.sliceOrder];
                    [temp[i], temp[i + 1]] = [temp[i + 1], temp[i]];
                    setImportState((s) => ({ ...s, sliceOrder: temp }));
                  }}
                >
                  <ArrowRight />
                </IconButton>
              </Grid>
            </Stack>
          );
        })}
      </Stack>
    </Stack>
  );
  const [height, setHeight] = useState(0);
  const [width, setWidth] = useState(0);

  useEffect(() => {
    if (!selectedImg) return;
    getImg(imgUrls[selectedImg]).then((img) => {
      setHeight(img.naturalHeight);
      setWidth(img.naturalWidth);
    });
  }, [imgUrls, selectedImg]);
  const [opacity, setOpacity] = useState(0.8);
  const alignments = selectedImg && currentImages[selectedImg]?.alignments;
  const handleXChange = (_: unknown, n: number) => {
    if (null === selectedImg) return;
    const newCurrImgs = [...currentImages];
    newCurrImgs[selectedImg].alignments.x = n;
    setCurrentImages(newCurrImgs);
  };
  const handleYChange = (_: unknown, n: number) => {
    if (null === selectedImg) return;
    const newCurrImgs = [...currentImages];
    newCurrImgs[selectedImg].alignments.y = n;
    setCurrentImages(newCurrImgs);
  };
  const handleRotZChange = (_: unknown, n: number) => {
    if (null === selectedImg) return;
    const newCurrImgs = [...currentImages];
    newCurrImgs[selectedImg].alignments.rotZ = n;
    setCurrentImages(newCurrImgs);
  };
  const mainImg = useRef<HTMLImageElement | null>(null);

  const naturalHeight = height;
  const naturalWidth = width;
  // console.log(naturalWidth, naturalHeight);

  const displayArea = alignments &&
    selectedImg &&
    selectedImg < currentImages.length && (
      <Grid
        style={{ width: "100%", position: "relative" }}
        direction="row"
        spacing={2}
        item
        container
        xs={12}
      >
        <Grid item xs={6}>
          <div
            style={{
              width: "500px",
              height: "500px",
              position: "relative",
              pointerEvents: "none",
            }}
          >
            <img
              src={imgUrls[selectedImg - 1]}
              style={{
                width: "100%",
                height: "100%",
                objectFit: "contain",
                position: "relative",
                transform: "",
              }}
            />
            <img
              ref={mainImg}
              src={imgUrls[selectedImg]}
              style={{
                position: "absolute",
                width: "100%",
                height: "100%",
                objectFit: "contain",
                left: 0,
                top: 0,
                opacity,
                transform:
                  "translateX(" +
                  (alignments.x / naturalWidth) * 100 +
                  "%)translateY(" +
                  (alignments.y / naturalHeight) * 100 +
                  "%)translateY(-50%)translateX(-50%)rotate(" +
                  alignments.rotZ +
                  "deg)translateY(50%)translateX(50%)" +
                  "",
              }}
            />
          </div>
        </Grid>

        <Grid item container sx={{ minWidth: 500 }} xs={6}>
          <Paper elevation={9} style={{ padding: 20 }}>
            <Box sx={{ width: 500 }}>
              <Typography id="input-slider" gutterBottom>
                Opacity
              </Typography>
              <Grid container spacing={2} alignItems="center">
                <Grid item>
                  <OpacityIcon />
                </Grid>
                <Grid item xs>
                  <Slider
                    value={opacity}
                    onChange={(_, n) => typeof n === "number" && setOpacity(n)}
                    aria-labelledby="input-slider"
                    step={0.01}
                    min={0}
                    max={1}
                  />
                </Grid>
                <Grid item>
                  <Input
                    style={{ width: 80 }}
                    value={opacity}
                    size="small"
                    onChange={(e) => setOpacity(parseFloat(e.target.value))}
                    inputProps={{
                      step: 10,
                      min: 0,
                      max: 100,
                      type: "number",
                      "aria-labelledby": "input-slider",
                    }}
                  />
                </Grid>
              </Grid>
            </Box>

            <Box sx={{ width: 500 }}>
              <Typography id="input-slider" gutterBottom>
                X Position
              </Typography>
              <Grid container spacing={2} alignItems="center">
                <Grid item>
                  <SwapHorizIcon />
                </Grid>
                <Grid item xs>
                  <Slider
                    onChange={(_, n) => handleXChange(null, n)}
                    aria-labelledby="input-slider"
                    min={-naturalWidth}
                    max={naturalWidth + Math.max(naturalWidth, naturalHeight)} //To handle all rotations
                    step={0.1}
                    value={alignments.x}
                  />
                </Grid>
                <Grid item>
                  <Input
                    style={{ width: 80 }}
                    value={alignments.x}
                    size="small"
                    onChange={(e) =>
                      handleXChange(null, parseFloat(e.target.value))
                    }
                    inputProps={{
                      step: 1,
                      min: -naturalWidth,
                      max: naturalWidth + Math.max(naturalWidth, naturalHeight),
                      type: "number",
                      "aria-labelledby": "input-slider",
                    }}
                  />
                </Grid>
              </Grid>
            </Box>

            <Box sx={{ width: 500 }}>
              <Typography id="input-slider" gutterBottom>
                Y Position
              </Typography>
              <Grid container spacing={2} alignItems="center">
                <Grid item>
                  <SwapVertIcon />
                </Grid>
                <Grid item xs>
                  <Slider
                    aria-labelledby="input-slider"
                    min={-naturalHeight - Math.max(naturalWidth, naturalHeight)} //To handle all rotation angles
                    max={naturalHeight}
                    step={0.1}
                    value={-alignments.y}
                    onChange={(_, n) => handleYChange(null, -n)}
                  />
                </Grid>
                <Grid item>
                  <Input
                    style={{ width: 80 }}
                    value={-alignments.y}
                    size="small"
                    onChange={(e) =>
                      handleYChange(null, -parseFloat(e.target.value))
                    }
                    inputProps={{
                      step: 1,
                      min: -naturalWidth,
                      max: naturalWidth + Math.max(naturalWidth, naturalHeight),
                      type: "number",
                      "aria-labelledby": "input-slider",
                    }}
                  />
                </Grid>
              </Grid>
            </Box>

            <Box sx={{ width: 500 }}>
              <Typography id="input-slider" gutterBottom>
                Rotation
              </Typography>
              <Grid container spacing={2} alignItems="center">
                <Grid item>
                  <LoopIcon />
                </Grid>
                <Grid item xs>
                  <Slider
                    aria-labelledby="input-slider"
                    min={-180}
                    max={180}
                    step={0.1}
                    value={alignments.rotZ}
                    onChange={handleRotZChange}
                  />
                </Grid>
                <Grid item>
                  <Input
                    style={{ width: 80 }}
                    value={alignments.rotZ}
                    size="small"
                    onChange={(e) =>
                      handleRotZChange(null, parseFloat(e.target.value))
                    }
                    inputProps={{
                      step: 0.1,
                      min: -180,
                      max: 180,
                      type: "number",
                      "aria-labelledby": "input-slider",
                    }}
                  />
                </Grid>
              </Grid>
            </Box>

            <Stack
              direction="row"
              spacing={2}
              style={{ width: "100%", paddingTop: 30 }}
            >
              <Button
                component="label"
                variant="contained"
                startIcon={<ClearIcon />}
                color="warning"
                onClick={(e) => {
                  handleXChange(null, 0);
                  handleYChange(null, 0);
                  handleRotZChange(null, 0);
                }}
              >
                Reset Alignment
              </Button>
              <div style={{ flex: "1 1 auto" }} />
              <Button component="label" variant="contained" startIcon={<ReplayIcon />}
                onClick={e => {
                  handleXChange(null, importState.pasteData[selectedImg]["px"] as number);
                  handleYChange(null, importState.pasteData[selectedImg]["py"] as number);
                  handleRotZChange(null, importState.pasteData[selectedImg]["theta"] as number);
                }}>
                Reload from Paste
              </Button>
            </Stack>
          </Paper>
        </Grid>

      </Grid>
    );

  return (
    <Stack style={{ padding: 20 }}>
      <Button variant="contained" component="label" startIcon={<ImageIcon />} style={{ textTransform: 'none' }}>
        Alignment Images ({numSlices} Required, {currentImages.length} Imported, {numSlices - currentImages.length} Missing)
        <input
          hidden
          multiple
          accept="image/*"
          type="file"
          onChange={(e) => {
            const files = e.target.files;
            if (!files) return;
            setCurrentImages(
              [...new Array(files.length)].map((_, i) => files[i]).map((f, i) => {
                  const alignments = importState.pasteData[i]
                    ? {
                        x: importState.pasteData[i]["px"] as number,
                        y: importState.pasteData[i]["py"] as number,
                        rotZ: importState.pasteData[i]["theta"] as number,
                      }
                    : { x: 0, y: 0, rotZ: 0 };
                  return {
                    file: f,
                    alignments: alignments,
                  };
                })
            );

            setImportState((s) => ({
              ...s,
              sliceOrder: [...new Array(numSlices)].map((_, i) => i),
            }));
          }}
        />
      </Button>
      {imageOrdering}
      {displayArea}
    </Stack>
  );
};
