/**
 * Classical Dirac Electron simulator (Barut–Zanghi model).
 * Port of CDE_4thRK.c — 4th-order Runge–Kutta time integration.
 * All constants are configurable via the params object.
 */

(function (global) {
  'use strict';

  function vec4(a, b, c, d) {
    return [a, b, c, d];
  }

  function copy4(src, dst) {
    dst[0] = src[0]; dst[1] = src[1]; dst[2] = src[2]; dst[3] = src[3];
    return dst;
  }

  function XDOT(zr, zi, xdot) {
    xdot[0] = zr[0]*zr[0] + zr[1]*zr[1] + zr[2]*zr[2] + zr[3]*zr[3] +
              zi[0]*zi[0] + zi[1]*zi[1] + zi[2]*zi[2] + zi[3]*zi[3];
    xdot[1] = 2 * (zr[0]*zr[3] + zi[0]*zi[3] + zr[1]*zr[2] + zi[1]*zi[2]);
    xdot[2] = 2 * (zr[0]*zi[3] - zi[0]*zr[3] - zr[1]*zi[2] + zi[1]*zr[2]);
    xdot[3] = 2 * (zr[0]*zr[2] + zi[0]*zi[2] - zr[1]*zr[3] - zi[1]*zi[3]);
  }

  function ZDOT(zr, zi, p, zrdot, zidot, lambda) {
    var l = lambda;
    zrdot[0] = (-p[0]*zi[0] + p[1]*zi[3] - p[2]*zr[3] + p[3]*zi[2]) / l;
    zrdot[1] = (-p[0]*zi[1] + p[1]*zi[2] + p[2]*zr[2] - p[3]*zi[3]) / l;
    zrdot[2] = ( p[0]*zi[2] - p[1]*zi[1] + p[2]*zr[1] - p[3]*zi[0]) / l;
    zrdot[3] = ( p[0]*zi[3] - p[1]*zi[0] - p[2]*zr[0] + p[3]*zi[1]) / l;
    zidot[0] = (-p[0]*zr[0] + p[1]*zr[3] + p[2]*zi[3] + p[3]*zr[2]) / (-l);
    zidot[1] = (-p[0]*zr[1] + p[1]*zr[2] - p[2]*zi[2] - p[3]*zr[3]) / (-l);
    zidot[2] = ( p[0]*zr[2] - p[1]*zr[1] - p[2]*zi[1] - p[3]*zr[0]) / (-l);
    zidot[3] = ( p[0]*zr[3] - p[1]*zr[0] + p[2]*zi[0] + p[3]*zr[1]) / (-l);
  }

  function PDOT(Fmunu, xdot, pdot, q) {
    pdot[0] =  q * (Fmunu[0][0]*xdot[0] + Fmunu[0][1]*xdot[1] + Fmunu[0][2]*xdot[2] + Fmunu[0][3]*xdot[3]);
    pdot[1] = -q * (Fmunu[1][0]*xdot[0] + Fmunu[1][1]*xdot[1] + Fmunu[1][2]*xdot[2] + Fmunu[1][3]*xdot[3]);
    pdot[2] = -q * (Fmunu[2][0]*xdot[0] + Fmunu[2][1]*xdot[1] + Fmunu[2][2]*xdot[2] + Fmunu[2][3]*xdot[3]);
    pdot[3] = -q * (Fmunu[3][0]*xdot[0] + Fmunu[3][1]*xdot[1] + Fmunu[3][2]*xdot[2] + Fmunu[3][3]*xdot[3]);
  }

  function ConstantEB(Fmunu, EX, EZ, BZ) {
    Fmunu[0][0] = 0;  Fmunu[0][1] = EX;  Fmunu[0][2] = 0;  Fmunu[0][3] = EZ;
    Fmunu[1][0] = -EX; Fmunu[1][1] = 0;  Fmunu[1][2] = -BZ; Fmunu[1][3] = 0;
    Fmunu[2][0] = 0;  Fmunu[2][1] = BZ;  Fmunu[2][2] = 0;  Fmunu[2][3] = 0;
    Fmunu[3][0] = -EZ; Fmunu[3][1] = 0;  Fmunu[3][2] = 0;  Fmunu[3][3] = 0;
  }

  function buildFmunu(EX, EZ, BZ) {
    var F = [
      [0, EX, 0, EZ],
      [-EX, 0, -BZ, 0],
      [0, BZ, 0, 0],
      [-EZ, 0, 0, 0]
    ];
    return F;
  }

  /**
   * Default parameters matching CDE_4thRK.c
   */
  function defaultParams() {
    return {
      lambda: 1,
      q: 1,
      QC: 0,
      EX: 0,
      EZ: 0,
      BZ: -0.9,
      dt: 1e-4,
      T: 30 * Math.PI,
      // initial momentum direction (radians): p_x = p0*cos(alpha), etc.
      initAlpha: 0.5 * Math.PI,
      initBeta: 0,
      initP0: 0  // 3-momentum scale; p0[0] then becomes sqrt(1 + p^2)
    };
  }

  /**
   * Initial state matching C code: x=0, zr=[1,0,0,-1], zi=0, p from alpha,beta
   */
  function initialState(params) {
    var p = params || defaultParams();
    var alpha = p.initAlpha, beta = p.initBeta, p0mag = p.initP0;
    var x0 = vec4(0, 0, 0, 0);
    var zr0 = vec4(1, 0, 0, -1);
    var zi0 = vec4(0, 0, 0, 0);
    var p0 = vec4(
      0,
      p0mag * Math.cos(alpha),
      p0mag * Math.sin(alpha) * Math.sin(beta),
      p0mag * Math.sin(alpha) * Math.cos(beta)
    );
    p0[0] = Math.sqrt(1 + p0[1]*p0[1] + p0[2]*p0[2] + p0[3]*p0[3]);
    return { x: x0, zr: zr0, zi: zi0, p: p0, t: 0 };
  }

  /**
   * Single RK4 step. Mutates state in place; returns the same state.
   */
  function step(state, params) {
    var lambda = params.lambda, q = params.q, dt = params.dt;
    var Fmunu = buildFmunu(params.EX, params.EZ, params.BZ);

    var x = state.x, zr = state.zr, zi = state.zi, p = state.p;
    var x1 = vec4(x[0], x[1], x[2], x[3]);
    var x2 = vec4(x[0], x[1], x[2], x[3]);
    var x3 = vec4(x[0], x[1], x[2], x[3]);
    var zr1 = vec4(zr[0], zr[1], zr[2], zr[3]);
    var zr2 = vec4(zr[0], zr[1], zr[2], zr[3]);
    var zr3 = vec4(zr[0], zr[1], zr[2], zr[3]);
    var zi1 = vec4(zi[0], zi[1], zi[2], zi[3]);
    var zi2 = vec4(zi[0], zi[1], zi[2], zi[3]);
    var zi3 = vec4(zi[0], zi[1], zi[2], zi[3]);
    var p1 = vec4(p[0], p[1], p[2], p[3]);
    var p2 = vec4(p[0], p[1], p[2], p[3]);
    var p3 = vec4(p[0], p[1], p[2], p[3]);

    var xdot = vec4(0,0,0,0), zrdot = vec4(0,0,0,0), zidot = vec4(0,0,0,0), pdot = vec4(0,0,0,0);
    var xdot1 = vec4(0,0,0,0), zrdot1 = vec4(0,0,0,0), zidot1 = vec4(0,0,0,0), pdot1 = vec4(0,0,0,0);
    var xdot2 = vec4(0,0,0,0), zrdot2 = vec4(0,0,0,0), zidot2 = vec4(0,0,0,0), pdot2 = vec4(0,0,0,0);
    var xdot3 = vec4(0,0,0,0), zrdot3 = vec4(0,0,0,0), zidot3 = vec4(0,0,0,0), pdot3 = vec4(0,0,0,0);

    XDOT(zr, zi, xdot);
    ZDOT(zr, zi, p, zrdot, zidot, lambda);
    PDOT(Fmunu, xdot, pdot, q);

    // k1 half step
    x1[0] += 0.5*dt*xdot[0];  x1[1] += 0.5*dt*xdot[1];  x1[2] += 0.5*dt*xdot[2];  x1[3] += 0.5*dt*xdot[3];
    zr1[0] += 0.5*dt*zrdot[0]; zr1[1] += 0.5*dt*zrdot[1]; zr1[2] += 0.5*dt*zrdot[2]; zr1[3] += 0.5*dt*zrdot[3];
    zi1[0] += 0.5*dt*zidot[0]; zi1[1] += 0.5*dt*zidot[1]; zi1[2] += 0.5*dt*zidot[2]; zi1[3] += 0.5*dt*zidot[3];
    p1[0] += 0.5*dt*pdot[0];  p1[1] += 0.5*dt*pdot[1];  p1[2] += 0.5*dt*pdot[2];  p1[3] += 0.5*dt*pdot[3];

    XDOT(zr1, zi1, xdot1);
    ZDOT(zr1, zi1, p1, zrdot1, zidot1, lambda);
    PDOT(Fmunu, xdot1, pdot1, q);

    x2[0] += 0.5*dt*xdot1[0]; x2[1] += 0.5*dt*xdot1[1]; x2[2] += 0.5*dt*xdot1[2]; x2[3] += 0.5*dt*xdot1[3];
    zr2[0] += 0.5*dt*zrdot1[0]; zr2[1] += 0.5*dt*zrdot1[1]; zr2[2] += 0.5*dt*zrdot1[2]; zr2[3] += 0.5*dt*zrdot1[3];
    zi2[0] += 0.5*dt*zidot1[0]; zi2[1] += 0.5*dt*zidot1[1]; zi2[2] += 0.5*dt*zidot1[2]; zi2[3] += 0.5*dt*zidot1[3];
    p2[0] += 0.5*dt*pdot1[0]; p2[1] += 0.5*dt*pdot1[1]; p2[2] += 0.5*dt*pdot1[2]; p2[3] += 0.5*dt*pdot1[3];

    XDOT(zr2, zi2, xdot2);
    ZDOT(zr2, zi2, p2, zrdot2, zidot2, lambda);
    PDOT(Fmunu, xdot2, pdot2, q);

    x3[0] += dt*xdot2[0]; x3[1] += dt*xdot2[1]; x3[2] += dt*xdot2[2]; x3[3] += dt*xdot2[3];
    zr3[0] += dt*zrdot2[0]; zr3[1] += dt*zrdot2[1]; zr3[2] += dt*zrdot2[2]; zr3[3] += dt*zrdot2[3];
    zi3[0] += dt*zidot2[0]; zi3[1] += dt*zidot2[1]; zi3[2] += dt*zidot2[2]; zi3[3] += dt*zidot2[3];
    p3[0] += dt*pdot2[0]; p3[1] += dt*pdot2[1]; p3[2] += dt*pdot2[2]; p3[3] += dt*pdot2[3];

    XDOT(zr3, zi3, xdot3);
    ZDOT(zr3, zi3, p3, zrdot3, zidot3, lambda);
    PDOT(Fmunu, xdot3, pdot3, q);

    x[0] += dt*(xdot[0] + 2*xdot1[0] + 2*xdot2[0] + xdot3[0])/6;
    x[1] += dt*(xdot[1] + 2*xdot1[1] + 2*xdot2[1] + xdot3[1])/6;
    x[2] += dt*(xdot[2] + 2*xdot1[2] + 2*xdot2[2] + xdot3[2])/6;
    x[3] += dt*(xdot[3] + 2*xdot1[3] + 2*xdot2[3] + xdot3[3])/6;
    zr[0] += dt*(zrdot[0] + 2*zrdot1[0] + 2*zrdot2[0] + zrdot3[0])/6;
    zr[1] += dt*(zrdot[1] + 2*zrdot1[1] + 2*zrdot2[1] + zrdot3[1])/6;
    zr[2] += dt*(zrdot[2] + 2*zrdot1[2] + 2*zrdot2[2] + zrdot3[2])/6;
    zr[3] += dt*(zrdot[3] + 2*zrdot1[3] + 2*zrdot2[3] + zrdot3[3])/6;
    zi[0] += dt*(zidot[0] + 2*zidot1[0] + 2*zidot2[0] + zidot3[0])/6;
    zi[1] += dt*(zidot[1] + 2*zidot1[1] + 2*zidot2[1] + zidot3[1])/6;
    zi[2] += dt*(zidot[2] + 2*zidot1[2] + 2*zidot2[2] + zidot3[2])/6;
    zi[3] += dt*(zidot[3] + 2*zidot1[3] + 2*zidot2[3] + zidot3[3])/6;
    p[0] += dt*(pdot[0] + 2*pdot1[0] + 2*pdot2[0] + pdot3[0])/6;
    p[1] += dt*(pdot[1] + 2*pdot1[1] + 2*pdot2[1] + pdot3[1])/6;
    p[2] += dt*(pdot[2] + 2*pdot1[2] + 2*pdot2[2] + pdot3[2])/6;
    p[3] += dt*(pdot[3] + 2*pdot1[3] + 2*pdot2[3] + pdot3[3])/6;

    state.t += dt;
    return state;
  }

  /**
   * Check if simulation should stop (time or mass shell violation).
   */
  function shouldStop(state, params) {
    if (state.t >= params.T) return true;
    var p = state.p;
    var m2 = p[0]*p[0] - p[1]*p[1] - p[2]*p[2] - p[3]*p[3];
    return m2 < -0.1;
  }

  /**
   * Spatial position (x, y, z) from state.
   */
  function position(state) {
    return { x: state.x[1], y: state.x[2], z: state.x[3] };
  }

  global.CDESim = {
    defaultParams: defaultParams,
    initialState: initialState,
    step: step,
    shouldStop: shouldStop,
    position: position
  };
})(typeof self !== 'undefined' ? self : this);
