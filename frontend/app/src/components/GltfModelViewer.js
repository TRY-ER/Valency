// import React, { Suspense, useRef } from 'react';
// import { Canvas } from '@react-three/fiber';
// import { OrbitControls, useGLTF, Environment } from '@react-three/drei';

// // This component loads and displays the GLTF model
// function Model({ modelPath, ...props }) {
//   const { scene } = useGLTF(modelPath);
//   // You can clone the scene if you want to reuse it multiple times
//   // or want to ensure modifications don't affect the original cache.
//   // For a single instance, directly using scene is fine.
//   return <primitive object={scene.clone()} {...props} />;
// }

// // Optional: A simple fallback loader
// function Loader() {
//   return (
//     <mesh>
//       <boxGeometry args={[1, 1, 1]} />
//       <meshStandardMaterial color="orange" />
//     </mesh>
//   );
//   // Or use Drei's Html loader:
//   // import { Html, useProgress } from '@react-three/drei';
//   // const { progress } = useProgress();
//   // return <Html center>{progress.toFixed(1)} % loaded</Html>;
// }

// function GltfModelViewer({ modelPath = '/models/DNA/dna.gltf' }) { // Default path
//   return (
//     <Canvas
//       camera={{ position: [0, 1.5, 5], fov: 50 }} // Adjust camera position and FOV as needed
//       shadows // Enable shadows in the scene
//       style={{ background: 'transparent' }} // Make canvas background transparent if needed
//     >
//       {/* Lights */}
//       <ambientLight intensity={0.7} />
//       <directionalLight
//         castShadow
//         position={[5, 10, 7.5]}
//         intensity={1.5}
//         shadow-mapSize-width={1024}
//         shadow-mapSize-height={1024}
//         shadow-camera-far={50}
//         shadow-camera-left={-10}
//         shadow-camera-right={10}
//         shadow-camera-top={10}
//         shadow-camera-bottom={-10}
//       />
//       <pointLight position={[-5, -5, -5]} intensity={0.5} />

//       {/* Environment for reflections and ambient lighting (optional but nice) */}
//       <Environment preset="sunset" background={false} /> {/* background={false} to keep canvas transparent */}

//       {/* Controls to interact with the model */}
//       <OrbitControls
//         enableDamping
//         dampingFactor={0.05}
//         autoRotate // Optional: make the model auto-rotate
//         autoRotateSpeed={0.5}
//       />

//       {/* Model Loading with Suspense */}
//       <Suspense fallback={<Loader />}>
//         <Model
//           modelPath={modelPath}
//           scale={[1, 1, 1]} // Adjust scale as needed
//           position={[0, 0, 0]} // Adjust position as needed
//           rotation={[0, 0, 0]} // Adjust rotation as needed
//         />
//       </Suspense>
//     </Canvas>
//   );
// }

// // Preload the model - good for performance if you know you'll use it
// useGLTF.preload('/models/DNA/dna.gltf'); // Adjust if you change the default path

// export default GltfModelViewer;


import React, { Suspense, useRef, useEffect } from 'react';
import { Canvas, useFrame } from '@react-three/fiber';
import { OrbitControls, useGLTF, Environment, Html, Grid } from '@react-three/drei';
import * as THREE from 'three'; // Import THREE for Box3 and Vector3

// This component loads and displays the GLTF model
function Model({ modelPath, ...props }) {
  const modelRef = useRef();
  // useGLTF loads the GLTF file.
  // It will automatically handle textures if they are defined in the GLTF's materials.
  const gltf = useGLTF(modelPath);

  // Log the loaded GLTF to inspect its structure
  useEffect(() => {
    if (gltf) {
      console.log("GLTF Loaded:", gltf);
      if (gltf.scene) {
        console.log("GLTF Scene:", gltf.scene);
        // Optional: Calculate and log bounding box to help with positioning/scaling
        const box = new THREE.Box3().setFromObject(gltf.scene);
        const size = box.getSize(new THREE.Vector3());
        const center = box.getCenter(new THREE.Vector3());
        console.log("Model Bounding Box Size:", size);
        console.log("Model Bounding Box Center:", center);

        // Example: If you want to auto-center and scale, you could do it here,
        // but it's often better to adjust the model in Blender or set scale/position props.
        // For instance, to roughly center:
        // if (modelRef.current) {
        //   modelRef.current.position.sub(center);
        // }
      }
    }
  }, [gltf]);

  // If the scene is not loaded, you might want to return null or a placeholder
  if (!gltf.scene) {
    console.warn(`GLTF scene not found for path: ${modelPath}`);
    return null;
  }

  // Use <primitive> to render the loaded GLTF scene
  // scene.clone() is good practice if you might reuse or modify the model instance
  return <primitive ref={modelRef} object={gltf.scene.clone()} {...props} />;
}

// A more visible fallback loader
function Loader() {
  return (
    <Html center>
      <div style={{ color: 'white', fontSize: '1.5em', textAlign: 'center' }}>
        Loading 3D Model...
      </div>
    </Html>
  );
}

function GltfModelViewer({ modelPath = '/models/DNA/dna.gltf' }) {
  return (
    <Canvas
      camera={{ position: [5, 5, 10], fov: 50, near: 0.1, far: 1000 }} // Adjusted camera
      shadows // Enable shadows in the scene
      style={{ background: 'transparent' }} // Keeps CSS background visible
      gl={{ antialias: true }} // Enable anti-aliasing
      onCreated={({ gl }) => {
        gl.toneMapping = THREE.ACESFilmicToneMapping; // Better color grading
        gl.outputEncoding = THREE.sRGBEncoding; // Correct color output
      }}
    >
      {/* Lights */}
      <ambientLight intensity={0.6} />
      <hemisphereLight skyColor={0xffffbb} groundColor={0x080820} intensity={0.5} />
      <directionalLight
        castShadow
        position={[10, 20, 15]} // Position the light to cast shadows effectively
        intensity={1.5}
        shadow-mapSize-width={2048} // Higher resolution shadows
        shadow-mapSize-height={2048}
        shadow-camera-far={50}
        shadow-camera-left={-15}
        shadow-camera-right={15}
        shadow-camera-top={15}
        shadow-camera-bottom={-15}
      />
      <pointLight position={[-10, -5, -10]} intensity={0.8} color="orange" />
      <pointLight position={[0, -5, 5]} intensity={0.3} />

      {/* Environment for reflections and overall scene lighting */}
      <Environment preset="city" background={false} />
      {/* Other presets: "sunset", "dawn", "night", "warehouse", "forest", "apartment", "studio", "city", "park", "lobby" */}

      {/* Controls to interact with the model */}
      <OrbitControls
        enableDamping
        dampingFactor={0.05}
        // autoRotate // Uncomment for auto-rotation
        // autoRotateSpeed={0.5}
        minDistance={1} // Prevent zooming in too close
        maxDistance={100} // Prevent zooming out too far
        target={[0, 1, 0]} // Point controls towards a typical model height
      />

      {/* Optional: Helper grid */}
      <Grid
        infiniteGrid
        cellSize={1}
        sectionSize={10}
        sectionColor={"#6f6f6f"}
        fadeDistance={50}
        fadeStrength={1}
      />

      {/* Ground plane to receive shadows */}
      <mesh rotation={[-Math.PI / 2, 0, 0]} position={[0, -0.01, 0]} receiveShadow>
        <planeGeometry args={[100, 100]} />
        <shadowMaterial opacity={0.3} />
        {/* You can use a standard material if you want a visible ground
        <meshStandardMaterial color="dimgray" /> */}
      </mesh>

      {/* Model Loading with Suspense */}
      <Suspense fallback={<Loader />}>
        <Model
          modelPath={modelPath}
          scale={[1, 1, 1]}       // START WITH 1,1,1 - ADJUST LATER!
          position={[0, 0, 0]}   // START WITH 0,0,0 - ADJUST LATER!
                                  // Models sometimes have an offset origin
          rotation={[0, 0, 0]}
          castShadow              // Ensure the model casts shadows
          receiveShadow           // Ensure the model can receive shadows (less common for main object)
        />
      </Suspense>
    </Canvas>
  );
}

// Preload the model - crucial for performance on first load
// Ensure this path matches the default or the one you intend to load first.
useGLTF.preload('/models/DNA/dna.gltf');

export default GltfModelViewer;