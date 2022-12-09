import {
  Button,
  Checkbox,
  FormControlLabel,
  FormGroup,
  Grid,
  Paper,
  Slider,
  Typography,
} from "@mui/material";
import { GizmoHelper, GizmoViewport, OrbitControls } from "@react-three/drei";
import { Canvas, useThree } from "@react-three/fiber";
import _ from "lodash";
import * as React from "react";
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import * as THREE from "three";
import { CurvesDisplay } from "./CurvesDisplay";
import "../api/threejsHeadeers";
import {
  datatype,
  pointToVector,
  colors as defaultColorArray,
  customFileExtension,
} from "../api/constants";
import { PointsDisplay } from "./PointsDisplay";
import { AreaDisplay } from "./AreaDisplay";
import { VolumeDisplay } from "./VolumeDisplay";
//Display List WebGL

const saveFileGeneric = (href: string, filename: string) => {
  // create "a" HTML element with href to file
  const link = document.createElement("a");
  link.href = href;
  link.download = filename;
  document.body.appendChild(link);
  link.click();

  // clean up "a" element & remove ObjectURL
  document.body.removeChild(link);
  URL.revokeObjectURL(href);
};

//Modified from https://stackoverflow.com/questions/55613438/reactwrite-to-json-file-or-export-download-no-server
const saveFile = (myData: Record<string, unknown>) => {
  // create file in browser
  const fileName = "savedata";
  const json = JSON.stringify(myData); //, null, 2);
  const blob = new Blob([json], { type: "application/json" });
  const href = URL.createObjectURL(blob);

  saveFileGeneric(href, fileName + customFileExtension);
};

const uid = Math.random().toString();

const Screengrabber = ({
  passthrough,
}: {
  passthrough: (takeScreenshot: () => string) => void;
}) => {
  const { gl, scene, camera } = useThree();

  const takeScreenshot = useCallback(() => {
    gl.render(scene, camera);
    return gl.domElement.toDataURL();
  }, [camera, gl, scene]);
  useEffect(() => passthrough(takeScreenshot), [takeScreenshot, passthrough]);
  return <mesh />;
};

export const CustomRenderer = ({
  data: data,
  setData,
  setIsOpen,
}: {
  data: datatype | undefined;
  setData: React.Dispatch<React.SetStateAction<datatype | undefined>>;
  setIsOpen: () => void;
}) => {
  const debounceRef = React.useRef<ReturnType<typeof setTimeout> | null>(null);
  //Data to display
  const [visuals, setVisuals] = useState({
    area: true,
    volume: true,
    contour: true,
    points: true,
  });

  const [doClusters, setDoClusters] = useState(false);

  const sliceNames = useMemo(
    () => (doClusters ? data?.sliceNames : data?.sliceNames),
    [data?.sliceNames, doClusters]
  );
  const slices = useMemo(() => (doClusters ? data?.slices : data?.slices), [
    data?.slices,
    doClusters,
  ]);
  const ctrs2Dvals = useMemo(
    () => (doClusters ? data?.ctrs2Dclusters : data?.ctrs2Dvals),
    [data?.ctrs2Dclusters, data?.ctrs2Dvals, doClusters]
  );
  const numFeatures = useMemo(
    () => (doClusters ? data?.nClusters : data?.nat),
    [data?.nClusters, data?.nat, doClusters]
  );
  const tris2Dvals = useMemo(
    () => (doClusters ? data?.tris2Dclusters : data?.tris2Dvals),
    [data?.tris2Dclusters, data?.tris2Dvals, doClusters]
  );
  const ctrs3Dvals = useMemo(
    () => (doClusters ? data?.ctrs3Dclusters : data?.ctrs3Dvals),
    [data?.ctrs3Dclusters, data?.ctrs3Dvals, doClusters]
  );
  const featureNames = useMemo(
    () =>
      doClusters
        ? [...new Array(data?.nClusters)].map((v, i) => "" + i)
        : data?.featureNames,
    [data?.featureNames, data?.nClusters, doClusters]
  );
  const values = useMemo(() => (doClusters ? data?.clusters : data?.values), [
    data?.clusters,
    data?.values,
    doClusters,
  ]);

  //Display User Settings
  const [activeSlices, setActiveSlices] = useState<
    { name: string; on: boolean }[]
  >([]);
  useEffect(() => {
    setActiveSlices(
      sliceNames?.map((v) => ({
        name: v,
        on: true,
      })) || []
    );
  }, [sliceNames]);

  const [activeGroups, setActiveGroups] = useState<
    { name: string; on: boolean }[]
  >([]);
  useEffect(() => {
    setActiveGroups(
      featureNames?.map((v, i) => ({
        name: i + 1 !== featureNames.length ? v : "No Tissue",
        on: i + 1 !== featureNames.length,
      })) || []
    );
  }, [featureNames]);

  const [colors, setColors] = useState(defaultColorArray);
  const [opacity, setOpacity] = useState(1);
  useEffect(() => {
    const num = activeGroups.length;
    setColors(defaultColorArray.slice(0, num));
  }, [activeGroups.length]);

  const camRef = useRef(null);
  const saveHandler = useCallback(() => {
    const d = {
      activeGroups,
      activeSlices,
      colors,
      data,
      visuals,
      opacity,
    };
    saveFile(d);
  }, [activeGroups, activeSlices, colors, data, opacity, visuals]);

  useEffect(() => {
    const stop = window.electronAPI.onSave(saveHandler);
    return () => {
      stop();
    };
  }, [saveHandler]);

  const openExisting = (
    <Button variant="contained" component="label" id={uid}>
      Open Geometry
      <input
        hidden
        multiple={false}
        accept={customFileExtension}
        type="file"
        onChange={async (e) => {
          const target = e.currentTarget as HTMLInputElement;
          const files = target?.files;
          const file = files?.[0];
          if (!file) return;
          const rawText = await file.text();
          const d = JSON.parse(rawText);
          setData(d.data || data);
          console.log(d);
          setTimeout(() => {
            //Placeholder to deal with sideffects of data change
            setActiveGroups(d.activeGroups || activeGroups);
            setActiveSlices(d.activeSlices || activeSlices);
            setColors(d.colors || colors);
            setVisuals(d.visuals || visuals);
            setOpacity(d.opacity !== undefined ? d.opacity : opacity);
          }, 1000);
        }}
      />
    </Button>
  );

  useEffect(() => {
    const openHandler = () => {
      document.getElementById(uid)?.click();
    };
    const stop = window.electronAPI.onOpen(openHandler);
    return () => {
      stop();
    };
  }, []);
  //Calculated Data
  const pointsBySlice = useMemo(
    () =>
      slices?.map((slice) =>
        slice.map((pt) => new THREE.Vector3(pt[0], pt[1], pt[2]))
      ),
    [slices]
  );

  const pointsData = useMemo(() => {
    const ptsByValue =
      pointsBySlice &&
      values
        ?.map(
          (slice) => slice.map((ptVals) => ptVals.indexOf(_.max(ptVals) || -1)) //Get all the max indices
        )
        .map((slice, sliceIndex) =>
          slice.map((primaryValue, ptIndex) => ({
            pt: pointsBySlice[sliceIndex][ptIndex],
            group: primaryValue,
          }))
        )
        .slice(1, -1)
        .filter((v, i) => activeSlices[i]?.on)
        .reduce((prev, current) => [...prev, ...current], []);
    return Object.entries(_.groupBy(ptsByValue, (v) => v.group)).map(
      (entry) => ({
        group: parseInt(entry[0]),
        pts: entry[1].map((v) => v.pt),
      })
    );
  }, [activeSlices, values, pointsBySlice]);

  const center = useMemo(() => {
    if (!pointsBySlice) return undefined;
    const totalNumberOfPoints = pointsBySlice.reduce(
      (p, slice) => p + slice.length,
      0
    );
    return pointsBySlice
      .reduce(
        (prev, slice) =>
          prev.add(
            slice.reduce((p, c) => p.add(c), new THREE.Vector3(0, 0, 0))
          ),
        new THREE.Vector3(0, 0, 0)
      )
      .divideScalar(totalNumberOfPoints || 1);
  }, [pointsBySlice]);

  const curvesFinalData = useMemo(() => {
    if (!ctrs2Dvals) return [];
    const ctrs = ctrs2Dvals
      .slice(1, -1)
      .filter((_, i) => activeSlices[i]?.on)
      .map((slice) =>
        slice
          .map((ctr) => ({
            points: ctr[0].map((p) => pointToVector(p)),
            segs: ctr[1],
          }))
          .map((ctr, i) => ({
            segments: ctr.segs.flatMap((v) => {
              return v.map((v1) => ctr.points[v1]);
            }),
            val: i,
          }))
      );

    return ctrs
      .reduce(
        (acc, slice) => {
          slice.forEach((ctr) => acc[ctr.val].push(...ctr.segments));
          return acc;
        },
        [...new Array(numFeatures)].map(() => [] as THREE.Vector3[])
      )
      .map((pts, group) => ({ pts, group }));
  }, [activeSlices, ctrs2Dvals, numFeatures]);
  const areaDisplayData = React.useMemo(() => {
    if (!tris2Dvals) return [];
    const ctrs = tris2Dvals
      .slice(1, -1)
      .filter((_, i) => activeSlices[i]?.on)
      .map((slice) => ({
        points: slice[0],
        tris: slice[1],
        vals: slice[2],
      }));

    const finalData = [...new Array(numFeatures)].map(
      () => [] as THREE.Vector3[]
    );

    ctrs.flatMap((slice) =>
      slice.tris
        .map((indices, indexIndex) => ({
          triangle: indices.map((i) => pointToVector(slice.points[i])),
          val: slice.vals[indexIndex],
        }))
        .forEach((triangle) =>
          finalData[triangle.val].push(...triangle.triangle)
        )
    );
    return finalData.map((pts, group) => ({ pts, group }));
  }, [activeSlices, numFeatures, tris2Dvals]);

  const volumes = useMemo(() => {
    function createPolygon(poly: THREE.Vector3[]) {
      //Assuming the polygon is a star
      const centerPoint = poly
        .reduce((a, b) => a.add(b), new THREE.Vector3(0, 0, 0))
        .divideScalar(poly.length);

      return poly.flatMap((e, i) => {
        return [poly[i], poly[(i + 1) % poly.length], centerPoint];
      });
    }
    if (!ctrs3Dvals) return [];
    return ctrs3Dvals
      .map((val) => ({
        points: val[0].map((p) => pointToVector(p)),
        rings: val[1],
      }))
      .map((val) =>
        val.rings.flatMap((indices) =>
          createPolygon(indices.map((index) => val.points[index]))
        )
      )
      .map((pts, group) => ({ pts, group }));
  }, [ctrs3Dvals]);
  const data2 = useMemo(
    () => ({
      volumes: volumes,
      areaDisplayData,
      curvesFinalData,
      pointsData,
    }),
    [areaDisplayData, curvesFinalData, pointsData, volumes]
  );
  const transformedData = useMemo(
    () => ({
      volumeDisplay: data2.volumes.filter((e) => activeGroups[e.group]?.on),
      areaDisplay: data2.areaDisplayData.filter(
        (e) => activeGroups[e.group]?.on
      ),
      curvesDisplay: data2.curvesFinalData.filter(
        (e) => activeGroups[e.group]?.on
      ),
      pointsDisplay: data2.pointsData.filter((e) => activeGroups[e.group]?.on),
    }),
    [
      activeGroups,
      data2.areaDisplayData,
      data2.curvesFinalData,
      data2.pointsData,
      data2.volumes,
    ]
  );

  const pointsDisplay = center && (
    <PointsDisplay
      center={center}
      groups={transformedData.pointsDisplay}
      colors={colors}
    />
  );

  const curvesDisplay = data && center && (
    <CurvesDisplay
      center={center}
      curvesFinalData={transformedData.curvesDisplay}
      colors={colors}
    />
  );

  const areaDisplay = data && center && (
    <AreaDisplay
      center={center}
      areaDisplayData={transformedData.areaDisplay}
      colors={colors}
      opacity={opacity}
    />
  );

  const volumeDisplay = data && center && (
    <VolumeDisplay
      center={center}
      volumes={transformedData.volumeDisplay}
      colors={colors}
      opacity={opacity}
    />
  );

  const canvas = useRef<(() => string) | null>(null);

  const saveCanvas = React.useCallback(() => {
    if (!canvas.current) return;
    const url = canvas.current();
    url && saveFileGeneric(url, "image");
  }, []);

  const renderSetup = (
    <>
      <Screengrabber passthrough={(v) => (canvas.current = v)} />
      <ambientLight />
      <OrbitControls makeDefault enableDamping={false} ref={camRef} />
      <GizmoHelper
        alignment="top-right" // widget alignment within scene
        margin={[80, 80]} // widget margins (X, Y)
      >
        <GizmoViewport
          axisColors={["red", "green", "blue"]}
          labelColor="white"
        />
      </GizmoHelper>
      <pointLight position={[10, 10, 10]} />
      <pointLight position={[-10, -10, -10]} />
    </>
  );

  const renderArea = (
    <Grid item xs={12}>
      <Canvas
        style={{
          height: 500,
          width: "100%",
          borderStyle: "solid",
          borderColor: "black",
          borderWidth: 3,
        }}
      >
        {renderSetup}

        {visuals.points && pointsDisplay}
        {visuals.contour && curvesDisplay}
        {visuals.area && areaDisplay}
        {visuals.volume && volumeDisplay}
      </Canvas>
    </Grid>
  );

  const primaryControlArea = (
    <Grid item container>
      <Grid item xs={3}>
        <FormGroup>
          <FormControlLabel
            control={<Checkbox checked={visuals.points} />}
            onChange={(_, checked) => {
              setVisuals((v) => ({ ...v, points: !!checked }));
            }}
            label={"Points"}
          />
        </FormGroup>
      </Grid>

      <Grid item xs={3}>
        <FormGroup>
          <FormControlLabel
            control={<Checkbox checked={visuals.contour} />}
            onChange={(_, checked) => {
              setVisuals((v) => ({ ...v, contour: !!checked }));
            }}
            label={"Contours"}
          />
        </FormGroup>
      </Grid>

      <Grid item xs={3}>
        <FormGroup>
          <FormControlLabel
            control={<Checkbox checked={visuals.area} />}
            onChange={(_, checked) => {
              setVisuals((v) => ({ ...v, area: !!checked }));
            }}
            label={"Areas"}
          />
        </FormGroup>
      </Grid>

      <Grid item xs={3}>
        <FormGroup>
          <FormControlLabel
            control={<Checkbox checked={visuals.volume} />}
            onChange={(_, checked) => {
              setVisuals((v) => ({ ...v, volume: !!checked }));
            }}
            label={"Volumes"}
          />
        </FormGroup>
      </Grid>
      <Grid item>
        <Typography>Volume Opacity</Typography>
        <Slider
          value={opacity}
          onChange={(_e, v) => setOpacity(v as number)}
          max={1}
          step={0.1}
          min={0}
        ></Slider>
      </Grid>
    </Grid>
  );

  const leftControlArea = (
    <Grid item container xs={12}>
      <Grid item>
        <Typography variant={"h5"}>Slices</Typography>
        <FormGroup>
          {activeSlices.map((active, i) => (
            <FormControlLabel
              key={i}
              control={<Checkbox checked={active.on} />}
              onChange={(_, checked) => {
                const oldData = activeSlices;
                oldData[i].on = checked;

                setActiveSlices([...oldData]);
              }}
              label={active.name}
            />
          ))}
        </FormGroup>
      </Grid>
    </Grid>
  );

  const rightControlArea = (
    <Grid item container xs={12}>
      <Grid item>
        <Typography variant={"h5"}>Groups</Typography>
        <FormGroup>
          {activeGroups.map((active, i) => (
            <div
              key={i}
              style={{
                display: "flex",
                alignItems: "center",
                width: "100%",
              }}
            >
              <div
                style={{
                  width: 80,
                  display: "flex",
                  alignItems: "center",
                  justifyContent: "space-around",
                  paddingRight: 10,
                }}
              >
                <Checkbox
                  checked={active.on}
                  onChange={(_, checked) => {
                    const oldData = activeGroups;
                    oldData[i].on = checked;
                    setActiveGroups([...oldData]);
                  }}
                />
                <input
                  style={{ width: 30 }}
                  type="color"
                  value={colors[i] as string}
                  onChange={(e) => {
                    const c = [...colors];
                    c[i] = e.target.value;
                    if (debounceRef.current) {
                      clearTimeout(debounceRef.current);
                    }
                    debounceRef.current = setTimeout(() => setColors(c), 100);
                  }}
                />
              </div>
              <div style={{ flex: 1 }}>
                <Typography>{active.name}</Typography>
              </div>
            </div>
          ))}
        </FormGroup>
      </Grid>
    </Grid>
  );

  return !data ? (
    <Grid container style={{ width: "100%" }}>
      <Button variant="outlined" onClick={() => setIsOpen()}>
        Import Data
      </Button>
      {openExisting}
    </Grid>
  ) : (
    <Grid container style={{ width: "100%" }}>
      <Button variant="outlined" onClick={() => setIsOpen()}>
        Import Data
      </Button>
      {openExisting}
      <Button onClick={saveHandler} variant="outlined">
        Save Geometry
      </Button>
      <Button onClick={() => setDoClusters((e) => !e)}>
        {doClusters ? "View Features" : "View Clusters"}
      </Button>
      <Button onClick={saveCanvas}>Download Picture</Button>
      <Grid item container xs={12} spacing={3}>
        <Grid item md={3} lg={2}>
          <Paper
            style={{ width: "100%", boxSizing: "border-box", padding: 15 }}
            elevation={4}
          >
            {leftControlArea}
          </Paper>
        </Grid>

        <Grid item md={3} lg={2}>
          <Paper
            style={{ width: "100%", boxSizing: "border-box", padding: 15 }}
            elevation={4}
          >
            {rightControlArea}
          </Paper>
        </Grid>
        <Grid item container md={6} lg={8}>
          <Paper
            style={{ width: "100%", boxSizing: "border-box", padding: 15 }}
            elevation={9}
          >
            {renderArea}
            {primaryControlArea}
          </Paper>
        </Grid>
      </Grid>
    </Grid>
  );
};
