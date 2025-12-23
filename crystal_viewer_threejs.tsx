import React, { useEffect, useRef, useState } from 'react';
import * as THREE from 'three';

const CrystalViewer = () => {
  const mountRef = useRef(null);
  const sceneRef = useRef(null);
  const rendererRef = useRef(null);
  const cameraRef = useRef(null);
  const animationRef = useRef(null);
  const [isRotating, setIsRotating] = useState(true);
  const [selectedStructure, setSelectedStructure] = useState('diamond');
  const [showBonds, setShowBonds] = useState(true);
  const [showUnitCell, setShowUnitCell] = useState(true);

  // Crystal structure data
  const crystalStructures = {
    diamond: {
      name: 'Diamond (C)',
      lattice: { a: 3.567, b: 3.567, c: 3.567, alpha: 90, beta: 90, gamma: 90 },
      atoms: [
        { element: 'C', position: [0, 0, 0], color: '#909090', radius: 0.7 },
        { element: 'C', position: [0.25, 0.25, 0.25], color: '#909090', radius: 0.7 },
        { element: 'C', position: [0.5, 0.5, 0], color: '#909090', radius: 0.7 },
        { element: 'C', position: [0.75, 0.75, 0.25], color: '#909090', radius: 0.7 },
        { element: 'C', position: [0.5, 0, 0.5], color: '#909090', radius: 0.7 },
        { element: 'C', position: [0.75, 0.25, 0.75], color: '#909090', radius: 0.7 },
        { element: 'C', position: [0, 0.5, 0.5], color: '#909090', radius: 0.7 },
        { element: 'C', position: [0.25, 0.75, 0.75], color: '#909090', radius: 0.7 },
      ]
    },
    perovskite: {
      name: 'Perovskite (CaTiO₃)',
      lattice: { a: 3.84, b: 3.84, c: 3.84, alpha: 90, beta: 90, gamma: 90 },
      atoms: [
        { element: 'Ca', position: [0, 0, 0], color: '#00ff00', radius: 1.0 },
        { element: 'Ti', position: [0.5, 0.5, 0.5], color: '#c0c0c0', radius: 0.6 },
        { element: 'O', position: [0.5, 0.5, 0], color: '#ff0000', radius: 0.6 },
        { element: 'O', position: [0.5, 0, 0.5], color: '#ff0000', radius: 0.6 },
        { element: 'O', position: [0, 0.5, 0.5], color: '#ff0000', radius: 0.6 },
      ]
    },
    nacl: {
      name: 'Rock Salt (NaCl)',
      lattice: { a: 5.64, b: 5.64, c: 5.64, alpha: 90, beta: 90, gamma: 90 },
      atoms: [
        { element: 'Na', position: [0, 0, 0], color: '#ab5cf2', radius: 1.02 },
        { element: 'Na', position: [0.5, 0.5, 0], color: '#ab5cf2', radius: 1.02 },
        { element: 'Na', position: [0.5, 0, 0.5], color: '#ab5cf2', radius: 1.02 },
        { element: 'Na', position: [0, 0.5, 0.5], color: '#ab5cf2', radius: 1.02 },
        { element: 'Cl', position: [0.5, 0, 0], color: '#1ff01f', radius: 1.81 },
        { element: 'Cl', position: [0, 0.5, 0], color: '#1ff01f', radius: 1.81 },
        { element: 'Cl', position: [0, 0, 0.5], color: '#1ff01f', radius: 1.81 },
        { element: 'Cl', position: [0.5, 0.5, 0.5], color: '#1ff01f', radius: 1.81 },
      ]
    },
    zincblende: {
      name: 'Zinc Blende (ZnS)',
      lattice: { a: 5.41, b: 5.41, c: 5.41, alpha: 90, beta: 90, gamma: 90 },
      atoms: [
        { element: 'Zn', position: [0, 0, 0], color: '#7d80b0', radius: 0.74 },
        { element: 'Zn', position: [0.5, 0.5, 0], color: '#7d80b0', radius: 0.74 },
        { element: 'Zn', position: [0.5, 0, 0.5], color: '#7d80b0', radius: 0.74 },
        { element: 'Zn', position: [0, 0.5, 0.5], color: '#7d80b0', radius: 0.74 },
        { element: 'S', position: [0.25, 0.25, 0.25], color: '#ffff30', radius: 1.04 },
        { element: 'S', position: [0.75, 0.75, 0.25], color: '#ffff30', radius: 1.04 },
        { element: 'S', position: [0.75, 0.25, 0.75], color: '#ffff30', radius: 1.04 },
        { element: 'S', position: [0.25, 0.75, 0.75], color: '#ffff30', radius: 1.04 },
      ]
    }
  };

  useEffect(() => {
    if (!mountRef.current) return;

    // Scene setup
    const scene = new THREE.Scene();
    scene.background = new THREE.Color(0x0a0a0a);
    sceneRef.current = scene;

    // Camera setup
    const camera = new THREE.PerspectiveCamera(
      50,
      mountRef.current.clientWidth / mountRef.current.clientHeight,
      0.1,
      1000
    );
    camera.position.set(15, 15, 15);
    camera.lookAt(0, 0, 0);
    cameraRef.current = camera;

    // Renderer setup
    const renderer = new THREE.WebGLRenderer({ antialias: true });
    renderer.setSize(mountRef.current.clientWidth, mountRef.current.clientHeight);
    renderer.setPixelRatio(window.devicePixelRatio);
    mountRef.current.appendChild(renderer.domElement);
    rendererRef.current = renderer;

    // Lighting
    const ambientLight = new THREE.AmbientLight(0xffffff, 0.5);
    scene.add(ambientLight);

    const directionalLight1 = new THREE.DirectionalLight(0xffffff, 0.8);
    directionalLight1.position.set(5, 5, 5);
    scene.add(directionalLight1);

    const directionalLight2 = new THREE.DirectionalLight(0xffffff, 0.4);
    directionalLight2.position.set(-5, -5, -5);
    scene.add(directionalLight2);

    // Mouse interaction
    let isDragging = false;
    let previousMousePosition = { x: 0, y: 0 };

    const onMouseDown = (e) => {
      isDragging = true;
      setIsRotating(false);
    };

    const onMouseUp = () => {
      isDragging = false;
    };

    const onMouseMove = (e) => {
      if (isDragging) {
        const deltaX = e.clientX - previousMousePosition.x;
        const deltaY = e.clientY - previousMousePosition.y;

        scene.rotation.y += deltaX * 0.01;
        scene.rotation.x += deltaY * 0.01;
      }
      previousMousePosition = { x: e.clientX, y: e.clientY };
    };

    const onWheel = (e) => {
      e.preventDefault();
      const delta = e.deltaY * 0.01;
      camera.position.multiplyScalar(1 + delta * 0.1);
    };

    renderer.domElement.addEventListener('mousedown', onMouseDown);
    renderer.domElement.addEventListener('mouseup', onMouseUp);
    renderer.domElement.addEventListener('mousemove', onMouseMove);
    renderer.domElement.addEventListener('wheel', onWheel);

    // Animation loop
    const animate = () => {
      animationRef.current = requestAnimationFrame(animate);
      
      if (isRotating) {
        scene.rotation.y += 0.005;
      }
      
      renderer.render(scene, camera);
    };
    animate();

    // Handle resize
    const handleResize = () => {
      if (!mountRef.current) return;
      camera.aspect = mountRef.current.clientWidth / mountRef.current.clientHeight;
      camera.updateProjectionMatrix();
      renderer.setSize(mountRef.current.clientWidth, mountRef.current.clientHeight);
    };
    window.addEventListener('resize', handleResize);

    // Cleanup
    return () => {
      window.removeEventListener('resize', handleResize);
      renderer.domElement.removeEventListener('mousedown', onMouseDown);
      renderer.domElement.removeEventListener('mouseup', onMouseUp);
      renderer.domElement.removeEventListener('mousemove', onMouseMove);
      renderer.domElement.removeEventListener('wheel', onWheel);
      if (animationRef.current) {
        cancelAnimationFrame(animationRef.current);
      }
      if (mountRef.current && renderer.domElement) {
        mountRef.current.removeChild(renderer.domElement);
      }
      renderer.dispose();
    };
  }, []);

  // Update structure when selection changes
  useEffect(() => {
    if (!sceneRef.current) return;

    const scene = sceneRef.current;
    
    // Clear existing structure
    while(scene.children.length > 0) { 
      const obj = scene.children[0];
      scene.remove(obj);
      if (obj.geometry) obj.geometry.dispose();
      if (obj.material) obj.material.dispose();
    }

    // Re-add lights
    const ambientLight = new THREE.AmbientLight(0xffffff, 0.5);
    scene.add(ambientLight);

    const directionalLight1 = new THREE.DirectionalLight(0xffffff, 0.8);
    directionalLight1.position.set(5, 5, 5);
    scene.add(directionalLight1);

    const directionalLight2 = new THREE.DirectionalLight(0xffffff, 0.4);
    directionalLight2.position.set(-5, -5, -5);
    scene.add(directionalLight2);

    const structure = crystalStructures[selectedStructure];
    const lattice = structure.lattice;

    // Create atoms
    structure.atoms.forEach((atom) => {
      const geometry = new THREE.SphereGeometry(atom.radius * 0.3, 32, 32);
      const material = new THREE.MeshPhongMaterial({
        color: atom.color,
        shininess: 80,
        specular: 0x444444
      });
      const sphere = new THREE.Mesh(geometry, material);
      
      const pos = atom.position;
      sphere.position.set(
        pos[0] * lattice.a - lattice.a / 2,
        pos[1] * lattice.b - lattice.b / 2,
        pos[2] * lattice.c - lattice.c / 2
      );
      
      scene.add(sphere);
    });

    // Create bonds
    if (showBonds) {
      const bondRadius = 0.08;
      const bondColor = 0x666666;
      
      for (let i = 0; i < structure.atoms.length; i++) {
        for (let j = i + 1; j < structure.atoms.length; j++) {
          const atom1 = structure.atoms[i];
          const atom2 = structure.atoms[j];
          
          const pos1 = new THREE.Vector3(
            atom1.position[0] * lattice.a - lattice.a / 2,
            atom1.position[1] * lattice.b - lattice.b / 2,
            atom1.position[2] * lattice.c - lattice.c / 2
          );
          
          const pos2 = new THREE.Vector3(
            atom2.position[0] * lattice.a - lattice.a / 2,
            atom2.position[1] * lattice.b - lattice.b / 2,
            atom2.position[2] * lattice.c - lattice.c / 2
          );
          
          const distance = pos1.distanceTo(pos2);
          
          // Only create bond if atoms are close enough
          if (distance < 2.5) {
            const direction = new THREE.Vector3().subVectors(pos2, pos1);
            const bondLength = direction.length();
            
            const geometry = new THREE.CylinderGeometry(bondRadius, bondRadius, bondLength, 8);
            const material = new THREE.MeshPhongMaterial({ color: bondColor });
            const bond = new THREE.Mesh(geometry, material);
            
            bond.position.copy(pos1.clone().add(direction.multiplyScalar(0.5)));
            bond.quaternion.setFromUnitVectors(
              new THREE.Vector3(0, 1, 0),
              direction.clone().normalize()
            );
            
            scene.add(bond);
          }
        }
      }
    }

    // Create unit cell
    if (showUnitCell) {
      const edges = [
        [[0, 0, 0], [1, 0, 0]],
        [[0, 0, 0], [0, 1, 0]],
        [[0, 0, 0], [0, 0, 1]],
        [[1, 0, 0], [1, 1, 0]],
        [[1, 0, 0], [1, 0, 1]],
        [[0, 1, 0], [1, 1, 0]],
        [[0, 1, 0], [0, 1, 1]],
        [[0, 0, 1], [1, 0, 1]],
        [[0, 0, 1], [0, 1, 1]],
        [[1, 1, 0], [1, 1, 1]],
        [[1, 0, 1], [1, 1, 1]],
        [[0, 1, 1], [1, 1, 1]],
      ];

      edges.forEach(([start, end]) => {
        const points = [
          new THREE.Vector3(
            start[0] * lattice.a - lattice.a / 2,
            start[1] * lattice.b - lattice.b / 2,
            start[2] * lattice.c - lattice.c / 2
          ),
          new THREE.Vector3(
            end[0] * lattice.a - lattice.a / 2,
            end[1] * lattice.b - lattice.b / 2,
            end[2] * lattice.c - lattice.c / 2
          )
        ];

        const geometry = new THREE.BufferGeometry().setFromPoints(points);
        const material = new THREE.LineBasicMaterial({ color: 0x00ffff, linewidth: 2 });
        const line = new THREE.Line(geometry, material);
        scene.add(line);
      });
    }

  }, [selectedStructure, showBonds, showUnitCell]);

  return (
    <div className="w-full h-screen bg-gray-950 flex flex-col">
      {/* Header */}
      <div className="bg-gradient-to-r from-blue-900 to-purple-900 p-4 shadow-lg">
        <h1 className="text-2xl font-bold text-white mb-2">Crystal Structure Viewer</h1>
        <p className="text-blue-200 text-sm">Professional Three.js-based crystallography visualization</p>
      </div>

      {/* Main content */}
      <div className="flex-1 flex">
        {/* Sidebar */}
        <div className="w-64 bg-gray-900 p-4 overflow-y-auto border-r border-gray-800">
          <div className="space-y-6">
            {/* Structure Selection */}
            <div>
              <h3 className="text-white font-semibold mb-3 text-sm uppercase tracking-wide">Structure</h3>
              <div className="space-y-2">
                {Object.entries(crystalStructures).map(([key, struct]) => (
                  <button
                    key={key}
                    onClick={() => setSelectedStructure(key)}
                    className={`w-full text-left px-3 py-2 rounded transition-colors ${
                      selectedStructure === key
                        ? 'bg-blue-600 text-white'
                        : 'bg-gray-800 text-gray-300 hover:bg-gray-700'
                    }`}
                  >
                    {struct.name}
                  </button>
                ))}
              </div>
            </div>

            {/* Display Options */}
            <div>
              <h3 className="text-white font-semibold mb-3 text-sm uppercase tracking-wide">Display</h3>
              <div className="space-y-2">
                <label className="flex items-center space-x-2 text-gray-300 cursor-pointer">
                  <input
                    type="checkbox"
                    checked={showBonds}
                    onChange={(e) => setShowBonds(e.target.checked)}
                    className="w-4 h-4"
                  />
                  <span>Show Bonds</span>
                </label>
                <label className="flex items-center space-x-2 text-gray-300 cursor-pointer">
                  <input
                    type="checkbox"
                    checked={showUnitCell}
                    onChange={(e) => setShowUnitCell(e.target.checked)}
                    className="w-4 h-4"
                  />
                  <span>Show Unit Cell</span>
                </label>
                <label className="flex items-center space-x-2 text-gray-300 cursor-pointer">
                  <input
                    type="checkbox"
                    checked={isRotating}
                    onChange={(e) => setIsRotating(e.target.checked)}
                    className="w-4 h-4"
                  />
                  <span>Auto Rotate</span>
                </label>
              </div>
            </div>

            {/* Structure Info */}
            <div>
              <h3 className="text-white font-semibold mb-3 text-sm uppercase tracking-wide">Lattice Parameters</h3>
              <div className="bg-gray-800 rounded p-3 text-sm">
                <div className="text-gray-400 space-y-1">
                  <div>a = {crystalStructures[selectedStructure].lattice.a} Å</div>
                  <div>b = {crystalStructures[selectedStructure].lattice.b} Å</div>
                  <div>c = {crystalStructures[selectedStructure].lattice.c} Å</div>
                  <div className="pt-2 border-t border-gray-700 mt-2">
                    <div>α = {crystalStructures[selectedStructure].lattice.alpha}°</div>
                    <div>β = {crystalStructures[selectedStructure].lattice.beta}°</div>
                    <div>γ = {crystalStructures[selectedStructure].lattice.gamma}°</div>
                  </div>
                </div>
              </div>
            </div>

            {/* Controls */}
            <div>
              <h3 className="text-white font-semibold mb-3 text-sm uppercase tracking-wide">Controls</h3>
              <div className="bg-gray-800 rounded p-3 text-xs text-gray-400 space-y-1">
                <div>🖱️ Drag to rotate</div>
                <div>🔄 Scroll to zoom</div>
                <div>⏸️ Click to pause rotation</div>
              </div>
            </div>
          </div>
        </div>

        {/* 3D Viewer */}
        <div className="flex-1 relative">
          <div ref={mountRef} className="w-full h-full" />
          
          {/* Overlay info */}
          <div className="absolute top-4 right-4 bg-black bg-opacity-70 text-white px-4 py-2 rounded-lg text-sm">
            <div className="font-semibold">{crystalStructures[selectedStructure].name}</div>
            <div className="text-gray-300 text-xs mt-1">
              {crystalStructures[selectedStructure].atoms.length} atoms
            </div>
          </div>
        </div>
      </div>

      {/* Footer */}
      <div className="bg-gray-900 border-t border-gray-800 px-4 py-2 text-xs text-gray-500 text-center">
        Powered by Three.js r128 • Materials Project-inspired architecture
      </div>
    </div>
  );
};

export default CrystalViewer;