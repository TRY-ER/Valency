import React, { useRef, useEffect, Suspense } from 'react';
import { Canvas, useFrame, useThree, useLoader } from '@react-three/fiber';
import { GLTFLoader } from 'three/examples/jsm/loaders/GLTFLoader';
import { OrbitControls, Environment, Loader, Stage, SpotLight, useGLTF } from '@react-three/drei';
import * as THREE from 'three';

// Preload the DNA model
useGLTF.preload('/models/DNA/dna.gltf');

function DnaModel() {
  // Use useGLTF instead of useLoader for better compatibility and error handling
  const { scene } = useGLTF('/models/DNA/dna.gltf');
  const modelRef = useRef();
  
  // Slow rotation effect
  useFrame(() => {
    if (modelRef.current) {
      modelRef.current.rotation.y += 0.001;
    }
  });

  // Scale and position the model
  useEffect(() => {
    if (modelRef.current) {
      modelRef.current.scale.set(0.6, 0.6, 0.6);
      modelRef.current.position.set(0, -2, 0); // Slightly lower position
      
      // Make sure model receives and casts shadows
      modelRef.current.traverse((child) => {
        if (child.isMesh) {
          child.castShadow = true;
          child.receiveShadow = true;
          child.material.needsUpdate = true;
          // Enhance material properties for better visibility
          if (child.material) {
            child.material.metalness = 0.3;
            child.material.roughness = 0.2;
          }
        }
      });
    }
  }, []);

  return (
    <primitive 
      ref={modelRef} 
      object={scene} 
      position={[0, 0, 0]} 
      castShadow
      receiveShadow
    />
  );
}

// A simple loading component that will be shown while the model loads
function LoadingScreen() {
  return (
    <div style={{ 
      position: 'absolute', 
      top: 0, 
      left: 0, 
      width: '100%', 
      height: '100%', 
      display: 'flex', 
      justifyContent: 'center', 
      alignItems: 'center',
      backgroundColor: 'rgba(0,0,0,0.1)',
      zIndex: 2
    }}>
      <div>Loading DNA model...</div>
    </div>
  );
}

function DnaBackground() {
  return (
    <div className="dna-background">
      <Canvas
        shadows
        camera={{ position: [5, 2, 8], fov: 45 }}
        style={{
          position: 'absolute',
          top: 0,
          left: 0,
          width: '100%',
          height: '100%',
          pointerEvents: 'none',
          zIndex: 1,
        }}
        gl={{ alpha: true, antialias: true }}
      >
        <color attach="background" args={['transparent']} />
        
        {/* Enhanced lighting setup */}
        <ambientLight intensity={0.8} />
        <directionalLight 
          position={[10, 10, 5]} 
          intensity={2}
          castShadow
          shadow-mapSize-width={1024}
          shadow-mapSize-height={1024}
        />
        <spotLight 
          position={[-10, 10, 5]} 
          intensity={1.5} 
          angle={0.6}
          penumbra={0.5}
          castShadow
        />
        <pointLight position={[0, 0, 5]} intensity={0.8} />
        
        {/* Environment adds natural lighting from an HDRI */}
        <Environment preset="sunset" />
        
        {/* Wrap the model in Suspense for loading handling */}
        <Suspense fallback={null}>
          <DnaModel />
        </Suspense>
        
        {/* Optional OrbitControls - useful for debugging but typically disabled in production */}
        {/* <OrbitControls
          enableZoom={false}
          enablePan={false}
          enableRotate={true}
          autoRotate={true}
          autoRotateSpeed={0.5}
        /> */}
      </Canvas>
      
      {/* External loader component from drei */}
      <Loader />
    </div>
  );
}

export default DnaBackground;
