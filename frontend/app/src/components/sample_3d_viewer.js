import React, { Suspense } from "react";
import ReactDOM from "react-dom";
import { Canvas } from "@react-three/fiber";
import { OrbitControls, Stats } from "@react-three/drei";

import { useGLTF } from "@react-three/drei";


function Plane(props) {
    // const { nodes, materials } = useGLTF("/models/small-airplane-v3.gltf");
    const { nodes, materials } = useGLTF("/models/DNA/dna.gltf");

    return (
        <group {...props}>
            <group scale={[0.5, 0.5, 0.5]}>
                <mesh
                    material={materials.White}
                    geometry={nodes["buffer-0-mesh-0"].geometry}
                />
                <mesh
                    material={materials.Red}
                    geometry={nodes["buffer-0-mesh-0_1"].geometry}
                />
                <mesh
                    material={materials.Gray}
                    geometry={nodes["buffer-0-mesh-0_2"].geometry}
                />
                <mesh
                    material={materials.Black}
                    geometry={nodes["buffer-0-mesh-0_3"].geometry}
                />
            </group>
        </group>
    );
}

useGLTF.preload("/models/DNA/dna.gltf");

const CanvasSample = () => {
    return (<>
        <Canvas style={{ height: 400, width: 800 }}>
            <pointLight position={[5, 5, 5]} />
            <Suspense fallback={null}>
                <Plane rotation={[0, Math.PI * 1.25, 0]} />
            </Suspense>
            <OrbitControls />
            <Stats />
        </Canvas>
    </>)
}

export default CanvasSample;