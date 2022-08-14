import React, { useEffect, useState } from "react";
import {
  Stack,
  Button,
  Grid,
  IconButton,
  Tooltip,
  Slider,
  Input,
  Typography,
} from "@mui/material";
import { ArrowLeft, ArrowRight } from "@mui/icons-material";
export const AlignmentPage = (): JSX.Element => {
  const [currentImages, setCurrentImages] = useState<
    { file: File; alignments: { x: number; y: number; rotZ: number } }[]
  >([]);
  const [selectedImg, setSelectedImg] = useState<number | null>(null);
  const [imgUrls, setImgUrls] = useState<string[]>([]);
  useEffect(() => {
    const func = async () => {
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
                  ? "solid 1px black"
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
  );
  const [opacity, setOpacity] = useState(0.5);
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
  const displayArea = alignments &&
    selectedImg &&
    selectedImg < currentImages.length && (
      <Grid
        style={{ width: "100%", position: "relative" }}
        direction="row"
        spacing={2}
      >
        <Grid item>
          <div
            style={{
              width: "100%",
              height: 700,
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
              }}
            />
            <img
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
                  alignments.x +
                  "%)translateY(" +
                  alignments.y +
                  "%)rotate(" +
                  alignments.rotZ +
                  "deg)",
              }}
            />
          </div>
        </Grid>
        <Grid item container sx={{ minWidth: 300 }} xs={12}>
          <Grid container item spacing={2} alignItems="center" xs={12}>
            <Typography>Opacity: </Typography>
            <Grid item xs>
              <Slider
                value={opacity}
                step={0.01}
                min={0}
                max={1}
                onChange={(_, n) => typeof n === "number" && setOpacity(n)}
              />
            </Grid>
          </Grid>

          {/* X Position */}
          <Grid container item spacing={2} alignItems="center" xs={12}>
            <Typography>X: </Typography>

            <Grid item xs>
              <Slider
                min={-100}
                max={100}
                step={0.1}
                value={alignments.x}
                onChange={handleXChange}
              />
            </Grid>
            <Grid item>
              <Input
                size="small"
                inputProps={{
                  step: 0.01,
                  min: -100,
                  max: 100,
                  type: "number",
                }}
                onChange={(e) =>
                  handleXChange(null, parseFloat(e.target.value))
                }
                value={alignments.x}
              />
            </Grid>
          </Grid>

          {/* Y Position */}
          <Grid container item spacing={2} alignItems="center" xs={12}>
            <Typography>Y: </Typography>

            <Grid item xs>
              <Slider
                min={-100}
                max={100}
                step={0.1}
                value={alignments.y}
                onChange={handleYChange}
              />
            </Grid>
            <Grid item>
              <Input
                size="small"
                inputProps={{
                  step: 0.01,
                  min: -100,
                  max: 100,
                  type: "number",
                }}
                onChange={(e) =>
                  handleYChange(null, parseFloat(e.target.value))
                }
                value={alignments.y}
              />
            </Grid>
          </Grid>

          {/* Z Rotation */}
          <Grid container item spacing={2} alignItems="center" xs={12}>
            <Typography>Rotation: </Typography>

            <Grid item xs>
              <Slider
                min={0}
                max={360}
                step={0.1}
                value={alignments.rotZ}
                onChange={handleRotZChange}
              />
            </Grid>
            <Grid item>
              <Input
                size="small"
                inputProps={{
                  step: 0.01,
                  min: 0,
                  max: 360,
                  type: "number",
                }}
                onChange={(e) =>
                  handleRotZChange(null, parseFloat(e.target.value))
                }
                value={alignments.rotZ}
              />
            </Grid>
          </Grid>
        </Grid>
      </Grid>
    );
  return (
    <Stack style={{ padding: 20 }}>
      <Button variant="contained" component="label">
        Alignment Images
        <input
          hidden
          multiple
          accept="image/*"
          type="file"
          onChange={(e) => {
            const files = e.target.files;

            if (!files) return;
            setCurrentImages(
              [...new Array(files.length)]
                .map((_, i) => files[i])
                .map((f) => ({ file: f, alignments: { x: 0, y: 0, rotZ: 0 } }))
            );
          }}
        />
      </Button>
      {imageOrdering}
      {displayArea}
    </Stack>
  );
};
