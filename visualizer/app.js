/**
 * Interactive 3D visualizer for Classical Dirac Electron trajectory.
 * Uses Three.js; trajectory grows in real time as the simulator steps.
 */

(function () {
  'use strict';

  var scene, camera, renderer, trajectoryLine, electronMesh, points;
  var state = null;
  var params = null;
  var running = false;
  var rafId = null;
  var cameraAngle = { theta: Math.PI / 4, phi: Math.PI / 4, radius: 14 };
  var mouseDown = false, lastMouse = { x: 0, y: 0 };
  var recording = false;
  var mediaRecorder = null;
  var recordedChunks = [];

  function getParams() {
    var dtExp = parseFloat(document.getElementById('dt').value);
    var dt = Math.pow(10, dtExp);
    return {
      lambda: parseFloat(document.getElementById('lambda').value),
      q: parseFloat(document.getElementById('q').value),
      QC: 0,
      EX: parseFloat(document.getElementById('EX').value),
      EZ: parseFloat(document.getElementById('EZ').value),
      BZ: parseFloat(document.getElementById('BZ').value),
      dt: dt,
      T: parseFloat(document.getElementById('T').value),
      initAlpha: parseFloat(document.getElementById('initAlpha').value),
      initBeta: parseFloat(document.getElementById('initBeta').value),
      initP0: parseFloat(document.getElementById('initP0').value)
    };
  }

  function stepsPerFrame() {
    return parseInt(document.getElementById('stepsPerFrame').value, 10) || 1;
  }

  function tailLengthSteps() {
    var units = parseInt(document.getElementById('tailLength').value, 10) || 8;
    return Math.max(2, units * 1000);
  }

  function syncSliderDisplay(id, value, decimals) {
    var el = document.getElementById(id + 'Val');
    if (el) el.textContent = typeof value === 'number' ? value.toFixed(decimals || 2) : value;
  }

  function initPanel() {
    document.getElementById('dt').addEventListener('input', function () {
      syncSliderDisplay('dt', Math.pow(10, parseFloat(this.value)), 4);
    });
    document.getElementById('stepsPerFrame').addEventListener('input', function () {
      document.getElementById('stepsPerFrameVal').textContent = this.value;
    });
    document.getElementById('lambda').addEventListener('input', function () {
      syncSliderDisplay('lambda', parseFloat(this.value), 2);
    });
    document.getElementById('q').addEventListener('input', function () {
      syncSliderDisplay('q', parseFloat(this.value), 2);
    });
    document.getElementById('EX').addEventListener('input', function () {
      syncSliderDisplay('EX', parseFloat(this.value), 2);
    });
    document.getElementById('EZ').addEventListener('input', function () {
      syncSliderDisplay('EZ', parseFloat(this.value), 2);
    });
    document.getElementById('BZ').addEventListener('input', function () {
      syncSliderDisplay('BZ', parseFloat(this.value), 2);
    });
    document.getElementById('initAlpha').addEventListener('input', function () {
      syncSliderDisplay('initAlpha', parseFloat(this.value), 2);
    });
    document.getElementById('initBeta').addEventListener('input', function () {
      syncSliderDisplay('initBeta', parseFloat(this.value), 2);
    });
    document.getElementById('initP0').addEventListener('input', function () {
      syncSliderDisplay('initP0', parseFloat(this.value), 2);
    });
    document.getElementById('tailLength').addEventListener('input', function () {
      document.getElementById('tailLengthVal').textContent = this.value;
    });
    syncSliderDisplay('dt', Math.pow(10, parseFloat(document.getElementById('dt').value)), 4);
    document.getElementById('stepsPerFrameVal').textContent = document.getElementById('stepsPerFrame').value;
    syncSliderDisplay('lambda', 1, 2);
    syncSliderDisplay('q', 1, 2);
    syncSliderDisplay('EX', 0, 2);
    syncSliderDisplay('EZ', 0, 2);
    syncSliderDisplay('BZ', -0.9, 2);
    syncSliderDisplay('initAlpha', 1.57, 2);
    syncSliderDisplay('initBeta', 0, 2);
    syncSliderDisplay('initP0', 0, 2);
    document.getElementById('tailLengthVal').textContent = document.getElementById('tailLength').value;
  }

  function initThree() {
    var canvas = document.getElementById('canvas');
    var width = canvas.clientWidth;
    var height = canvas.clientHeight;

    scene = new THREE.Scene();
    scene.background = new THREE.Color(0x0c0c10);

    camera = new THREE.PerspectiveCamera(50, width / height, 0.1, 1000);
    camera.position.set(8, 8, 8);
    camera.lookAt(0, 0, 0);

    renderer = new THREE.WebGLRenderer({ canvas: canvas, antialias: true });
    renderer.setSize(width, height);
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));

    // Comet tail: last N points, fading from tail (alpha 0) to head (alpha 1)
    var tailGeom = new THREE.BufferGeometry();
    tailGeom.setAttribute('position', new THREE.Float32BufferAttribute([], 3));
    tailGeom.setAttribute('alpha', new THREE.Float32BufferAttribute([], 1));
    var tailMaterial = new THREE.ShaderMaterial({
      transparent: true,
      depthWrite: false,
      vertexShader: [
        'attribute float alpha;',
        'varying float vAlpha;',
        'void main() {',
        '  vAlpha = alpha;',
        '  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);',
        '}'
      ].join('\n'),
      fragmentShader: [
        'varying float vAlpha;',
        'void main() {',
        '  gl_FragColor = vec4(0.51, 0.55, 0.98, vAlpha);',
        '}'
      ].join('\n')
    });
    trajectoryLine = new THREE.Line(tailGeom, tailMaterial);
    scene.add(trajectoryLine);

    // Current electron position (sphere)
    var sphereGeom = new THREE.SphereGeometry(0.08, 16, 12);
    var sphereMat = new THREE.MeshBasicMaterial({ color: 0xfbbf24 });
    electronMesh = new THREE.Mesh(sphereGeom, sphereMat);
    scene.add(electronMesh);

    // Grid helper (optional, subtle)
    var grid = new THREE.GridHelper(20, 20, 0x27272a, 0x1f1f23);
    scene.add(grid);

    window.addEventListener('resize', function () {
      var w = canvas.clientWidth;
      var h = canvas.clientHeight;
      camera.aspect = w / h;
      camera.updateProjectionMatrix();
      renderer.setSize(w, h);
    });

    canvas.addEventListener('mousedown', function (e) {
      mouseDown = true;
      lastMouse.x = e.clientX;
      lastMouse.y = e.clientY;
    });
    window.addEventListener('mousemove', function (e) {
      if (!mouseDown) return;
      var dx = e.clientX - lastMouse.x;
      var dy = e.clientY - lastMouse.y;
      cameraAngle.theta -= dx * 0.01;
      cameraAngle.phi = Math.max(0.05, Math.min(Math.PI - 0.05, cameraAngle.phi + dy * 0.01));
      lastMouse.x = e.clientX;
      lastMouse.y = e.clientY;
    });
    window.addEventListener('mouseup', function () { mouseDown = false; });
    canvas.addEventListener('wheel', function (e) {
      e.preventDefault();
      cameraAngle.radius = Math.max(2, Math.min(80, cameraAngle.radius + e.deltaY * 0.05));
    }, { passive: false });
  }

  function updateCamera() {
    var r = cameraAngle.radius;
    var t = cameraAngle.theta;
    var p = cameraAngle.phi;
    camera.position.set(r * Math.sin(p) * Math.cos(t), r * Math.sin(p) * Math.sin(t), r * Math.cos(p));
    camera.lookAt(0, 0, 0);
    camera.updateProjectionMatrix();
  }

  function updateTailGeometry() {
    var L = tailLengthSteps();
    var n = Math.floor(points.length / 3);
    if (n < 2) {
      trajectoryLine.geometry.setDrawRange(0, 0);
      return;
    }
    var take = Math.min(L, n);
    var start = (n - take) * 3;
    var posArr = new Float32Array(take * 3);
    var alphaArr = new Float32Array(take);
    for (var i = 0; i < take; i++) {
      posArr[i * 3] = points[start + i * 3];
      posArr[i * 3 + 1] = points[start + i * 3 + 1];
      posArr[i * 3 + 2] = points[start + i * 3 + 2];
      alphaArr[i] = take > 1 ? i / (take - 1) : 1;
    }
    trajectoryLine.geometry.setAttribute('position', new THREE.BufferAttribute(posArr, 3));
    trajectoryLine.geometry.setAttribute('alpha', new THREE.BufferAttribute(alphaArr, 1));
    trajectoryLine.geometry.setDrawRange(0, take);
    trajectoryLine.geometry.attributes.position.needsUpdate = true;
    trajectoryLine.geometry.attributes.alpha.needsUpdate = true;
  }

  function resetSimulation() {
    params = getParams();
    state = CDESim.initialState(params);
    points = [];
    var pos = CDESim.position(state);
    points.push(pos.x, pos.y, pos.z);
    updateTailGeometry();
    electronMesh.position.set(pos.x, pos.y, pos.z);
    updateStatus(false);
  }

  function updateStatus(stopped) {
    var msg = stopped ? 'Stopped. Press Run.' : 'Running… t = ' + (state ? state.t.toFixed(3) : '0');
    var el = document.getElementById('status');
    if (el) el.textContent = msg;
  }

  function animate() {
    rafId = requestAnimationFrame(animate);
    updateCamera();
    if (state) params = getParams();
    if (!running || !state || !params) {
      if (points && points.length >= 6) updateTailGeometry();
      renderer.render(scene, camera);
      return;
    }
    var n = stepsPerFrame();
    for (var i = 0; i < n; i++) {
      if (CDESim.shouldStop(state, params)) {
        running = false;
        document.getElementById('playPause').textContent = 'Run';
        updateStatus(true);
        break;
      }
      CDESim.step(state, params);
      var pos = CDESim.position(state);
      points.push(pos.x, pos.y, pos.z);
      electronMesh.position.set(pos.x, pos.y, pos.z);
    }
    updateTailGeometry();
    if (state) updateStatus(false);
    renderer.render(scene, camera);
  }

  function playPause() {
    if (!state) resetSimulation();
    params = getParams();
    running = !running;
    document.getElementById('playPause').textContent = running ? 'Pause' : 'Run';
    updateStatus(!running);
  }

  function reset() {
    running = false;
    if (rafId) cancelAnimationFrame(rafId);
    document.getElementById('playPause').textContent = 'Run';
    resetSimulation();
    rafId = requestAnimationFrame(animate);
  }

  function startRecording() {
    if (!window.MediaRecorder) {
      alert('Recording not supported in this browser. Try Chrome or Firefox.');
      return;
    }
    var canvas = renderer.domElement;
    var stream = canvas.captureStream(30);
    var mimeType = 'video/webm';
    if (MediaRecorder.isTypeSupported('video/webm;codecs=vp9')) {
      mimeType = 'video/webm;codecs=vp9';
    } else if (MediaRecorder.isTypeSupported('video/webm;codecs=vp8')) {
      mimeType = 'video/webm;codecs=vp8';
    }
    recordedChunks = [];
    var options = { videoBitsPerSecond: 2500000, mimeType: mimeType };
    try {
      mediaRecorder = new MediaRecorder(stream, options);
    } catch (e) {
      mediaRecorder = new MediaRecorder(stream);
    }
    mediaRecorder.ondataavailable = function (e) {
      if (e.data && e.data.size > 0) recordedChunks.push(e.data);
    };
    mediaRecorder.onstop = function () {
      var blob = new Blob(recordedChunks, { type: mimeType });
      var url = URL.createObjectURL(blob);
      var a = document.createElement('a');
      a.href = url;
      a.download = 'classical_dirac_electron_' + new Date().toISOString().slice(0, 19).replace(/:/g, '-') + '.webm';
      a.click();
      URL.revokeObjectURL(url);
      recordedChunks = [];
      recording = false;
      document.getElementById('recordBtn').textContent = 'Record';
    };
    mediaRecorder.start(1000);
    recording = true;
    document.getElementById('recordBtn').textContent = 'Stop recording';
  }

  function stopRecording() {
    if (mediaRecorder && recording && mediaRecorder.state !== 'inactive') {
      mediaRecorder.stop();
    }
  }

  document.getElementById('playPause').addEventListener('click', playPause);
  document.getElementById('reset').addEventListener('click', reset);
  document.getElementById('recordBtn').addEventListener('click', function () {
    if (recording) stopRecording();
    else startRecording();
  });

  initPanel();
  initThree();
  resetSimulation();
  rafId = requestAnimationFrame(animate);
})();
