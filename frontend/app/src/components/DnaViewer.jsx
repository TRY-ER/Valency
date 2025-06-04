// import React, { Suspense, useRef } from 'react';
// import { Canvas, useFrame } from '@react-three/fiber';
// import { useGLTF, OrbitControls, Environment } from '@react-three/drei';

// // Preload the DNA model
// useGLTF.preload('/models/DNA/dna.gltf');

// function DnaModel() {
//   const { scene } = useGLTF('/models/DNA/dna.gltf');
//   const modelRef = useRef();
  
//   // Rotate the model
//   useFrame(() => {
//     if (modelRef.current) {
//       modelRef.current.rotation.y += 0.005; // Slightly faster rotation for visibility
//     }
//   });

//   return (
//     <primitive 
//       ref={modelRef} 
//       object={scene} 
//       scale={[0.8, 0.8, 0.8]}
//       position={[0, 0, 0]} 
//       castShadow
//       receiveShadow
//     />
//   );
// }

// function DnaViewer() {
//   return (
//     <div className="dna-model-container">
//       <Canvas
//         shadows
//         camera={{ position: [0, 0, 5], fov: 75 }}
//         style={{
//           width: '100%',
//           height: '100%',
//           borderRadius: '12px',
//         }}
//       >
//         {/* Clear background */}
//         {/* <color attach="background" args={['#00000000']} /> */}
        
//         {/* Strong lighting setup for visibility */}
//         <ambientLight intensity={1.5} />
//         <spotLight position={[10, 10, 10]} angle={0.3} penumbra={1} intensity={2} castShadow />
//         <pointLight position={[-10, -10, -10]} intensity={1} />
        
//         {/* Environment for realistic lighting */}
//         <Environment preset="city" />
        
//         {/* The DNA model */}
//         <Suspense fallback={null}>
//           <DnaModel />
//         </Suspense>
        
//         {/* Controls for auto-rotation */}
//         <OrbitControls 
//           enableZoom={true} 
//           enablePan={true} 
//           enableRotate={true}
//           autoRotate={true}
//           autoRotateSpeed={1}
//         />
//       </Canvas>
//     </div>
//   );
// }

// export default DnaViewer;


