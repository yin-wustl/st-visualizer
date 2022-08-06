import { colors, defaultColor, shrink } from "../api/constants";
import * as THREE from "three";
import { Points } from "./threejsComponents/Points";
import * as React from "react";

export const PointsDisplay = ({
  center,
  groups,
}: {
  center: THREE.Vector3;
  groups: {
    group: number;
    pts: THREE.Vector3[];
  }[];
}): JSX.Element => {
  const translate = center
    ? ([-center.x, -center.y, -center.z] as const)
    : ([0, 0, 0] as const);

  return (
    <>
      {Object.values(groups).map((g) => {
        return (
          <Points
            key={g.group}
            points={g.pts}
            translate={translate}
            scale={shrink}
            color={colors[g.group] || defaultColor}
          />
        );
      })}
    </>
  );
};