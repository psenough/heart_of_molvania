#version 300 es

precision highp float;
precision highp int;
precision highp sampler2DArray;

uniform sampler2DArray t0;
uniform sampler2DArray t1;
uniform sampler2DArray t2;
uniform int layer;
uniform float u_timer;
uniform float u_delta_timer;
uniform vec2 u_resolution;

in vec2 v_st;

out vec4 color;
float iTime;
vec2 uv;
float PI = 3.14159265359;

/*float checkerBool (return float(h.x>.5==h.y>.5);}
float checkerBool2(vec2 b=gthv(.5)floatbool2
float checkerBoolT(vec2 b=gthv(cos(iTime)*.45+.5)floatbool2*/
float checkerSign(vec2 v){ v=sign(mod(v,2.)-.5+cos(iTime)*.4); return v.x*v.y; }

mat2 rot2(float r){float c=cos(r),s=sin(r);return mat2(c,s,-s,c);}

float checkerSin(vec2 v){return sign(cos(v.x)*sin(v.y));}



vec2 Rot(vec2 p, float t) {
    float c = cos(t); float s = sin(t);
    return vec2(p.x*c+p.y*s,
                -p.x*s+p.y*c);
}
vec2 RotCS(vec2 p, float c, float s) {
    return vec2( p.x*c+p.y*s,
                -p.x*s+p.y*c);
}

/* Created by Nikita Miropolskiy, nikat/2013
    * This work is licensed under a 
    * Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
    * http://creativecommons.org/licenses/by-nc-sa/3.0/
    *  - You must attribute the work in the source code 
    *    (link to https://www.shadertoy.com/view/XsX3zB).
    *  - You may not use this work for commercial purposes.
    *  - You may distribute a derivative work only under the same license.
    */

/* discontinuous pseudorandom uniformly distributed in [-0.5, +0.5]^3 */
vec3 random3(vec3 c) {
    float j = 4096.0*sin(dot(c,vec3(17.0, 59.4, 15.0)));
    vec3 r;
    r.z = fract(512.0*j);
    j *= .125;
    r.x = fract(512.0*j);
    j *= .125;
    r.y = fract(512.0*j);
    r = r-0.5;
    
    //rotate for extra flow!
    float t = -iTime*.5;
    r.xy = Rot(r.xy,t);

    return r;
}

/* skew constants for 3d simplex functions */
const float F3 =  0.3333333;
const float G3 =  0.1666667;

/* 3d simplex noise */
float noise(vec3 p) {
        /* 1. find current tetrahedron T and its four vertices */
        /* s, s+i1, s+i2, s+1.0 - absolute skewed (integer) coordinates of T vertices */
        /* x, x1, x2, x3 - unskewed coordinates of p relative to each of T vertices*/
        
        /* calculate s and x */
        vec3 s = floor(p + dot(p, vec3(F3)));
        vec3 x = p - s + dot(s, vec3(G3));
        
        /* calculate i1 and i2 */
        vec3 e = step(vec3(0.0), x - x.yzx);
        vec3 i1 = e*(1.0 - e.zxy);
        vec3 i2 = 1.0 - e.zxy*(1.0 - e);
        
        /* x1, x2, x3 */
        vec3 x1 = x - i1 + G3;
        vec3 x2 = x - i2 + 2.0*G3;
        vec3 x3 = x - 1.0 + 3.0*G3;
        
        /* 2. find four surflets and store them in d */
        vec4 w, d;
        
        /* calculate surflet weights */
        w.x = dot(x, x);
        w.y = dot(x1, x1);
        w.z = dot(x2, x2);
        w.w = dot(x3, x3);
        
        /* w fades from 0.6 at the center of the surflet to 0.0 at the margin */
        w = max(0.6 - w, 0.0);
        
        /* calculate surflet components */
        d.x = dot(random3(s), x);
        d.y = dot(random3(s + i1), x1);
        d.z = dot(random3(s + i2), x2);
        d.w = dot(random3(s + 1.0), x3);
        
        /* multiply d by w^4 */
        w *= w;
        w *= w;
        d *= w;
        
        /* 3. return the sum of the four surflets */
        return dot(d, vec4(52.0));
}


//iq 2d simplex noise

vec2 hash( vec2 p )
{
    p = vec2( dot(p,vec2(127.1,311.7)),
                dot(p,vec2(269.5,183.3)) );

    vec2 h = -1.0 + 2.0*fract(sin(p)*43758.5453123);

#if 1	
    //extra rotations for more flow!
    float t = -iTime*0.7;
    float co = cos(t); float si = sin(t);	
    h = RotCS(h,co,si);
#endif
    return h;
}


float noise( in vec2 p )
{
    const float K1 = 0.366025404; // (sqrt(3)-1)/2;
    const float K2 = 0.211324865; // (3-sqrt(3))/6;

    vec2 i = floor( p + (p.x+p.y)*K1 );
    
    vec2 a = p - i + (i.x+i.y)*K2;
    vec2 o = (a.x>a.y) ? vec2(1.0,0.0) : vec2(0.0,1.0); //vec2 of = 0.5 + 0.5*vec2(sign(a.x-a.y), sign(a.y-a.x));
    vec2 b = a - o + K2;
    vec2 c = a - 1.0 + 2.0*K2;

#if 1	
    //even more extra rotations for more flow!
    float t = iTime*.5;
    float co = cos(t); float si = sin(t);	
    a = RotCS(a,co,si);
    b = RotCS(b,co,si);
    c = RotCS(c,co,si);
#endif
    
    vec3 h = max( 0.5-vec3(dot(a,a), dot(b,b), dot(c,c) ), 0.0 );

    vec3 n = h*h*h*h*vec3( dot(a,hash(i+0.0)), dot(b,hash(i+o)), dot(c,hash(i+1.0)));

    return dot( n, vec3(70.0) );
    
}

float pot(vec2 pos)
{
    float t = iTime*.1;

    vec3 p = vec3(pos+vec2(iTime*.4,0.),t);
    
    float n = noise(p);
    n += 0.5 *noise(p*2.13);
    n += 3. * noise(pos*0.333);
    
    return n;
}

vec2 field(vec2 pos)
{
    float s = 1.5;
    pos *= s;
    
    float n = pot(pos);
    
    float e = 0.1;
    float nx = pot(vec2(pos+vec2(e,0.)));
    float ny = pot(vec2(pos+vec2(0.,e)));
    
    return vec2(-(ny-n),nx-n)/e;
}





vec2 iResolution = vec2(640.,480.);
int STEPS = 30;
float FAR = 60.0;
float PIXELR = 320.;
int BOUNCES = 3;
float SAMPLES = 4.0;
float CTIME = 0.0;

//Hash methods from https://www.shadertoy.com/view/4djSRW
//#define HASHSCALE3 vec3(.1031, .1030, .0973)
vec3 HASHSCALE3 = vec3(443.897, 441.423, 437.195);
vec3 hash33(vec3 p3){
    p3 = fract(p3 * HASHSCALE3);
    p3 += dot(p3, p3.yxz+19.19);
    return fract((p3.xxy + p3.yxx)*p3.zyx);
}

float HASHSCALE1 = 443.8975;
float hash13(vec3 p3){
    p3  = fract(p3 * HASHSCALE1);
    p3 += dot(p3, p3.yzx + 19.19);
    return fract((p3.x + p3.y) * p3.z);
}


//Distance functions from Mercury's SDF library
//http://mercury.sexy/hg_sdf/

// Maximum/minumum elements of a vector
float vmax3(vec3 v) {
    return max(max(v.x, v.y), v.z);
}

// Box: correct distance to corners
float fBox(vec3 p, vec3 b) {
    vec3 d = abs(p) - b;
    return length(max(d, vec3(0))) + vmax3(min(d, vec3(0)));
}

// Rotate around a coordinate axis (i.e. in a plane perpendicular to that axis) by angle <a>.
// Read like this: R(p.xz, a) rotates "x towards z".
// This is fast if <a> is a compile-time constant and slower (but still practical) if not.
void pR(inout vec2 p, float a) {
    p = cos(a)*p + sin(a)*vec2(p.y, -p.x);
}

// 3D noise function (IQ)
float noise2(vec3 p){
    vec3 ip = floor(p);
    p -= ip;
    vec3 s = vec3(7.0,157.0,113.0);
    vec4 h = vec4(0.0, s.yz, s.y+s.z)+dot(ip, s);
    p = p*p*(3.0-2.0*p);
    h = mix(fract(sin(h)*43758.5), fract(sin(h+s.x)*43758.5), p.x);
    h.xy = mix(h.xz, h.yw, p.y);
    return mix(h.x, h.y, p.z);
}

vec2 dist(vec3 p){
    vec3 pp = p;
    pR(pp.xy, CTIME);
    float tunnel = -fBox(pp, vec3(4.0, 4.0, 2.0*FAR));
    
    pp = p;
    pR(pp.xz, CTIME*0.75);
    pR(pp.yz, CTIME*0.25);
    float box = fBox(pp, vec3(0.5))-noise2((pp+CTIME*0.5)*0.5)*2.0;
    
    float scene = min(tunnel, box);
    float id = 0.0;
    
    if(box < tunnel){
        id = 1.0;
    }
    
    return vec2(scene, id);
}

vec3 normals(vec3 p){
    vec3 eps = vec3(PIXELR, 0.0, 0.0);
    return normalize(vec3(
        dist(p+eps.xyy).x-dist(p-eps.xyy).x,
        dist(p+eps.yxy).x-dist(p-eps.yxy).x,
        dist(p+eps.yyx).x-dist(p-eps.yyx).x
    ));
}

//Enhanced sphere tracing algorithm introduced by Mercury

// Sign function that doesn't return 0
float sgn(float x) {
    return (x < 0.0)?-1.0:1.0;
}

vec2 march(vec3 ro, vec3 rd){
    float t = 0.001;//EPSILON;
    float step = 0.0;

    float omega = 1.0;//muista testata eri arvoilla! [1,2]
    float prev_radius = 0.0;

    float candidate_t = t;
    float candidate_error = 1000.0;
    float sg = sgn(dist(ro).x);

    vec3 p = vec3(0.0);

    for(int i = 0; i < STEPS; ++i){
        p = rd*t+ro;
        float sg_radius = sg*dist(p).x;
        float radius = abs(sg_radius);
        step = sg_radius;
        bool fail = omega > 1. && (radius+prev_radius) < step;
        if(fail){
            step -= omega * step;
            omega = 1.;
        }
        else{
            step = sg_radius*omega;
        }
        prev_radius = radius;
        float error = radius/t;

        if(!fail && error < candidate_error){
            candidate_t = t;
            candidate_error = error;
        }

        if(!fail && error < PIXELR || t > FAR){
            break;
        }
        t += step;
    }
    //discontinuity reduction
    float er = candidate_error;
    for(int j = 0; j < 6; ++j){
        float radius = abs(sg*dist(p).x);
        p += rd*(radius-er);
        t = length(p-ro);
        er = radius/t;

        if(er < candidate_error){
            candidate_t = t;
            candidate_error = er;
        }
    }
    if(t <= FAR || candidate_error <= PIXELR){
        t = candidate_t;
    }
    
    p = ro+rd*t;
    float id = dist(p).y;
    
    return vec2(t, id);
}

//returns material of the object hit
// emissive color is xyz, and reflectance w
vec4 getMaterial(float obj, vec3 p){
    vec3 base = vec3(0.1);
    float reflectance = 0.0;
    float m = mod(p.z-(CTIME*10.0), 8.0) - 4.0;
    vec3 col = vec3(0.5, 0.4, 0.8);
    
    if(obj == 0.0){
        if(m > 0.0 && m > 2.0){
            base = col;
        }
        else if( m < 0.0 && m > -2.0){
            base = col.rgr;
        }
        reflectance = m > 0.0 ? 0.2 : 0.5;
    }
    else if(obj == 1.0){
        base = col.brg*0.2;
        reflectance = 0.4;
    }
    
    
    return vec4(base, reflectance);
}

vec3 render(vec3 o, vec3 d, vec2 uv){
    
    vec3 ro = o;
    vec3 rd = d;
    
    vec3 pixel_color = vec3(0.0);
    vec3 absorption_factor = vec3(1.0);
    
    for(int i = 0; i < BOUNCES; ++i){
        vec2 t = march(ro, rd);
        vec3 p = ro+rd*t.x;
        
        if(t.y < 0.0 || t.x > FAR){
            break;
        }
        
        //material.xyz == emissive
        //material.w == reflectance
        vec4 material = getMaterial(t.y, p);
        pixel_color += material.xyz * absorption_factor;
        absorption_factor *= material.w;
        
        vec3 n = normals(p);
        ro = p+(n*0.02);
        if(t.y == 0.0){
            rd = reflect(rd,n);
            //Thanks to fizzer to introducing this skew thing! :)
            rd = normalize(rd + (hash33(vec3(uv, float(i))) - 0.5)*0.1); 
        }
        else if(t.y == 1.0){
            rd = reflect(rd,n);
        }
        
    }
    
    return pixel_color;
}

// distance functions for the heart tunnel
// original shader from https://www.shadertoy.com/view/wsSBRc
float distCustom(float x, float y) {
    float n = 0.45 * abs(x) + y;
    return log(
        x * x + 1.8 * n * n
    );
}
float fn(float x, float multi, float offset) {
    return max(0.0, min(1.0, (sin(x * 3.14159265) + offset) * multi));
}
float fn2(float x) {
    return max(0.0, (sin(x) + 0.3) * 0.8);
}


float trange(float a, float b) {
    return max(0., min(1., (u_timer - a) / (b - a)));
}
float easeIn(float t) {
    return t * t * t;
}
float easeOut(float t) {
    float f = t - 1.0;
    return f * f * f + 1.0;
}
float sqr(float x) { 
    return x*x; 
}

void main() {
    vec4 color0 = vec4(0.);

    //int stuff = 0;
    
    //float sync = 2400.0;
    //if (u_timer - sync * floor(u_timer/sync) > sync*0.5) stuff = 1;
    
    float period = 475.;
    int panit = int(uint(floor(u_timer / period)) % uint(5));

    float fadeIn = min(1., u_timer / 10200.);
    float fadeOut = max(0., (u_timer - 120000.) / (133735.-120000.));
    
    float heartPulseRaw = cos((u_timer/1000.-31.825)/(0.860625) *3.14*2.)*.5 + .5;
    float heartPulse = pow(heartPulseRaw, 4.);

    float heartZoomRate = 0.00025;
    iTime = u_timer*heartZoomRate;	

    uv = vec2(v_st.xy-0.5)*2.0;
    uv.x *= u_resolution.x/u_resolution.y;
    //uv *= 0.5;
    //uv = vec2(uv.x*cos(iTime)-uv.y*sin(iTime), uv.x*sin(iTime)+uv.y*cos(iTime));

    float pinch = fadeOut;

    
    if (u_timer > 56000./*65000.*/ && u_timer < 78600.) {
        pinch = sqr(sin(trange(56000., 78600.)*PI));
    }
    if (u_timer > 84000. && u_timer < 105310.) {
        float tt = trange(84000., 105310.);
        pinch = heartPulseRaw*heartPulseRaw * sqr(sin(tt*PI)) * 0.5;
    }

    if (pinch > 0.) {
        float r = length(uv);
        float a = atan(uv.y, uv.x);
        r = mix(r, sqrt(r)*1.0, pinch);
        uv = r * vec2(cos(a), sin(a));
    }

    // sine distort
    if (u_timer < 10200.) {
        float t = easeOut(u_timer / 10200.);
        uv.x += sin(uv.y*PI + u_timer/1000.) *0.25*(1.-t);
        uv.y += cos(uv.x*PI + u_timer/1250.+PI/2.4663) *0.125*(1.-t);
    }

    if (u_timer > 75000. && u_timer < 78600.) { // horz row distort
        float tween = (u_timer - 75000.)/(78600.-75000.);
        float ramp = pow( sin( tween * PI), 2.);
        float rows = mix(32., 256., tween);//32.+32.*ramp;
        float amp = 0.125 * ramp;
        float freq = 8.;//2. + ramp * 8.;
        float row = floor(v_st.y*rows) / rows;
        //vec2 baseDistortion = vec2( hash22(vec2(row+params.phase, 0.0)).x-0.5, 0.0 ) * params.amplitude;
        uv.x += hash13(vec3(row, 0., 0.)) /*sin(row*PI*freq)*/ * amp;
        //uv = rot2(tween*PI*2.*tween*4.) * uv;
    }
    else if (u_timer > 77600.) { //} && u_timer < 103000.) { // vert row distort
        //float tween = (u_timer - 78600.)/(103000.-78600.);
        //float ramp = pow( sin( tween * PI), 2.);
        float rows = 256.;//32.+32.*ramp;
        float amp = max(0., sin(u_timer*0.00785))*0.075;// 0.025;//0.5 * ramp;
        float freq = 8.;//2. + ramp * 8.;
        float row = floor(v_st.x*rows) / rows;
        //vec2 baseDistortion = vec2( hash22(vec2(row+params.phase, 0.0)).x-0.5, 0.0 ) * params.amplitude;
        vec2 hh = hash(vec2(u_timer/1000000., row));
        if (hh.x > 0.9)
            uv.y += hh.y * amp;//hash13(vec3(row, 0., 0.)) /*sin(row*PI*freq)*/ * amp;
    }


    const float PI = 3.14159265;
    const float density = 1.5;
    float shift = sin(iTime)*0.15 + 0.15;
    vec3 color1 = vec3(1.0, shift, shift);
    vec3 color2 = vec3(1.0, 0.5, shift);

    uv *= 2.0;

    if (u_timer > 105000.) {
        float mx = (u_timer - 105000.)/10000.;
        mx *= mx; 
        mx *= mx;
        uv *= 1. + mx;
        iTime += mix(0., 4096.*heartZoomRate, mx);//mx*heartZoomRate;
    }
    
    float dist = log(uv.x*uv.x+uv.y*uv.y) / 2.;
    float distH = distCustom(uv.x, uv.y);
    float angle = atan(uv.y, uv.x);
    
    float timeH = -iTime;

    const float T_BRIDGE = 26622.; // end of beating heart, time slow down
    const float T_BRIDGE_END = 31700.;

    // slowed down sync hack section
    if ((u_timer > T_BRIDGE) && (u_timer < T_BRIDGE_END)) {
        float tween = (u_timer - T_BRIDGE) / (T_BRIDGE_END - T_BRIDGE);
        float dt = u_timer - T_BRIDGE;
        iTime = T_BRIDGE*heartZoomRate + dt*mix(heartZoomRate, heartZoomRate/8., tween);
        timeH = -iTime;
    }
    
    // heart tunnel, original shadertoy from HaleyHalcyon
    vec2 muv = vec2( cos(0.5-uv.x + sin(uv.x*uv.y + iTime)), sin(0.5+uv.y) + cos(0.5+uv.y) + 0.15*sin(sin(uv.x*20.0)*0.2));
    vec3 col2 = floor(mod(muv.xxx + (sin(iTime*0.1)-muv.yxx)*0.8-iTime*0.25,0.15)*50.0);
    //col2 = mix(vec3(0.0), col2, fadeIn);

    // Time varying pixel color

    float offset1 = 0.1, offset2 = 1.1;
    if ((u_timer > T_BRIDGE_END) && (u_timer < 76400.)) {
        offset2 += heartPulse*.2;
    }

    float c1 = fn((distH + timeH) * density + offset1, 1.1, -0.3) * col2.y; // heart rings 1
    float c2 = fn((distH + timeH) * density + offset2, 1.5, -0.3) * col2.x; // heart rings 2
    //float c3 = fn(dist * 4.0 + angle / PI + iTime * speed + PI, 3.0, -0.8);
    
    if (u_timer > 76000.) {
        float t = trange(76000., 90600.);
        float tt = trange(90000., 120000.);
        c1 = mix(c1, fn2(distH * 3.5 + angle + timeH * 4./*speed*/ + PI)*mix(col2.y,1.,1.-tt), t);
        c2 = mix(c2, fn2(distH * 4.5 - angle + timeH * 4./*speed*/ * -1.5 + PI)*mix(col2.x, 1., tt*0.25), t);
    }
    

    // twirl, original shadertoy from Ivan Weston
    float r = length(uv);			
    float sum = 0.0;
    iTime = u_timer*0.00005;
    for (int i = 0 ; i < 16; i++) {
        if (i < 16+int(sin(iTime)*16.0)) {
            float theta1 = (5.0*atan(uv.y, uv.x)-r*PI*4.0*cos(float(i)+iTime))+ cos(iTime);
            float awesome = pow(clamp(1.0-acos(cos(theta1)), 0.0, 1.0), PI);
            sum += awesome;
        }
    }
    vec4 twirl = vec4(1.);			
    twirl.r = cos(sum*1.0+cos(iTime*1.0)*2.0)*.5+.5;
    twirl.g = cos(sum*1.0+cos(iTime*2.0))*.5+.5;
    twirl.b = cos(sum*1.0+cos(iTime*3.0))*.5+.5;

/*		
    // plasma
    vec4 plasma = vec4(0.5+0.5*sin(iTime + vec3(1.0+sin(iTime + uv.y),2,uv.y+0.5+0.5*sin(uv.x*20.0+5.0*cos(uv.x+uv.y+iTime)))),1.);

    // swirl, original shadertoy from Antonalog
    vec2 old_uv = uv;
    //float lod = 0.;			
    vec3 d = vec3(0.);
    vec3 e = vec3(0.);
    for (int i=0; i<25; i++)
    {
        //texture(t2, vec3(v_st, i));
        d += texture(t0,vec3(uv+u_timer*0.00005,0)).xyz;
        e += texture(t0,vec3(-uv.yx*3.+u_timer*0.000125,0)).xyz;
        
        vec2 new_uv = field(uv)*.00625*.5;
    
        //lod += length(new_uv)*5.;
        uv += new_uv;
    }
    vec3 c = texture(t0,vec3(uv*.1+u_timer*0.000025,0)).xyz;
    d *= (1./50.);
    e *= (1./50.);
    c = mix(c,d,length(d));
    c = mix(c,e,length(e));
    vec4 swirl = vec4(c,1.);

    
    
    // rimina's tunnel blob thing
    //vec2 uv = fragCoord.xy / iResolution.xy;
    vec2 q = old_uv; //-1.0+2.0*uv;
    //q.x *= iResolution.x/iResolution.y;
    q *= period*0.01;
    vec3 ro = vec3(0.0, 0.0, 5.0);
    vec3 rt = vec3(0.0, 0.0, -2.0);
    vec3 z = normalize(rt-ro);
    vec3 x = normalize(cross(z, vec3(0.0, 1.0, 0.0)));
    vec3 y = normalize(cross(x, z));
    vec3 color1 = vec3(0.0);
    for(float i = 0.0; i < SAMPLES; ++i){
        //from http://www.iquilezles.org/www/articles/simplepathtracing/simplepathtracing.htm
        CTIME = u_timer*0.0005 + 0.6*(1.0/24.0)*hash13(vec3(uv, u_timer*0.0005));
        
        vec3 rd = normalize(mat3(x, y, z)*vec3(q, radians(60.0)));
        color1 += render(ro, rd+sin(iTime)*0.25, uv);
    }
    color1 /= SAMPLES;
    color1 = smoothstep(0.2, 0.9, color1);
    color1 = pow(color1, 1.0/vec3(2.2));
    vec4 rimina = vec4(color1, 1.0);


    // get mixed colors of pancake_1 and distance to twirl
    vec4 mixedv1 = color0;
    vec4 col1[5];
    float dist1[5];
    for (int i=0; i < dist1.length(); i++) {
        //if (stuff == 0) col[i] = texture(t1, vec3(v_st, i));
        //	else col[i] = texture(t0, vec3(v_st, i));
        col1[i] = texture(t2, vec3(v_st, i));
        dist1[i] = distance(twirl, col1[i]);
        //mixedv = mix(mixedv, col[i], sin(iTime*100.0+float(i+panit))*0.5+0.5);
        mixedv1 = mix(mixedv1, col1[i], abs(sin(iTime*5.+float(i)))*0.5);
    }
    
    // get mixed colors of pancake_2 and distance to twirl
    vec4 mixedv2 = color0;
    vec4 col2[4];
    float dist2[4];
    for (int i=0; i < dist2.length(); i++) {
        //if (stuff == 0) col[i] = texture(t1, vec3(v_st, i));
        //	else col[i] = texture(t0, vec3(v_st, i));
        col2[i] = texture(t0, vec3(v_st, i));
        dist2[i] = distance(twirl, col2[i]);
        //mixedv = mix(mixedv, col[i], sin(iTime*100.0+float(i+panit))*0.5+0.5);
        mixedv2 = mix(mixedv2, col2[i], abs(sin(iTime*8.+float(i)))*2.0);
    }
    
    // get mixed colors of tpolm and distance to plasma
    vec4 mixedv3 = color0;
    vec4 col3[11];
    float dist3[11];
    for (int i=0; i < dist3.length(); i++) {
        //if (stuff == 0) col[i] = texture(t1, vec3(v_st, i));
        //	else col[i] = texture(t0, vec3(v_st, i));
        col3[i] = texture(t1, vec3(v_st, i));
        dist3[i] = distance(plasma, col3[i]);
        //mixedv = mix(mixedv, col[i], sin(iTime*100.0+float(i+panit))*0.5+0.5);
        mixedv3 = mix(mixedv3, col3[i], abs(sin(iTime*5.+float(i)))*2.5);
    }

    // calculate shortest distance (in color) from pancake_2 to twirl
    float shortest_dist = 0.0;
    float second_shortest_dist = 0.0;
    int shortest_idx = 1;
    vec4 shortest_idx_v = vec4(0.0);
    vec4 second_shortest_idx_v = vec4(0.0);
    float td = distance(twirl, col1[0]);
    shortest_dist = td;
    second_shortest_dist = td;
    for (int i=0; i < dist1.length(); i++) {
        td = dist1[i];
        if (td <= shortest_dist) {
            shortest_idx = i;
            second_shortest_dist = shortest_dist;
            second_shortest_idx_v = col1[i];
            shortest_idx_v = col1[i];
        } else if (td <= second_shortest_dist) {
            second_shortest_dist = td;
            second_shortest_idx_v = col1[i];
        }
    }
    
    // calculate shortest distance (in color) from tpolm to plasma
    float shortest_dist2 = 0.0;
    float second_shortest_dist2 = 0.0;
    int shortest_idx2 = 1;
    vec4 shortest_idx_v2 = vec4(0.0);
    vec4 second_shortest_idx_v2 = vec4(0.0);
    float td2 = distance(swirl, col3[0]);
    shortest_dist2 = td2;
    second_shortest_dist2 = td2;
    for (int i=0; i < dist3.length(); i++) {
        td2 = dist3[i];
        if (td2 <= shortest_dist2) {
            shortest_idx2 = i;
            second_shortest_dist2 = shortest_dist2;
            second_shortest_idx_v2 = col3[i];
            shortest_idx_v2 = col3[i];
        } else if (td2 <= second_shortest_dist2) {
            second_shortest_dist2 = td2;
            second_shortest_idx_v2 = col3[i];
        }
    }

    // checkerboard splitting twirl and pancake_1, original shadertoy from ollj
    vec4 checkerboard = vec4(1.);
    float t = iTime*.01;
    vec2 v=uv;
    v*=rot2(t*0.1);
    v+=rot2(t*9.)[0]*32.;//position.xy //f.x+=cos(t*9.)*32.;//position.x //f.y+=sin(t*9.)*32.;//position.y
    v*=cos(t*.01)*.5+.51;//scale.xy
    v*=8.0;//scale.xy
    float s = checkerSign(v);
    s = s/0.4+.5;
    if (s>.5) checkerboard = mixedv1;
        else checkerboard = col1[uint(shortest_idx+panit) % uint(5)];


    // timeline
    if (u_timer < 900.) {
        color0 = mixedv1;
    } else if ((u_timer > 900.) && (u_timer < 4200.)) {
        color0 = texture(t2, vec3(v_st, panit));
        color0 *= 1.0 - (mod(u_timer, period) / period);
    } else if ((u_timer > 4400.) && (u_timer < 7850.)) {
        color0 = texture(t2, vec3(v_st, panit));
        color0 *= 1.0 - (mod(u_timer, period) / period);
    } else if ((u_timer > 8100.) && (u_timer < 21200.)) {
        color0 = twirl;
        if ( distance( vec4(0.,0.,0.,1.0), color0 ) > 0.3) {
            color0 = shortest_idx_v; //texture(t2, vec3(v_st, panit));
            if ( (u_timer < 11500.) ||
                    ((u_timer > 11750.) && (u_timer < 15200.)) ||
                    ((u_timer > 15650.) && (u_timer < 18800.)) ||
                    (u_timer > 19350.)
                ) color0 *= 1.5 - (mod(u_timer, period) / period);
        }
    } else if ((u_timer > 21200.) && (u_timer < 52000.)) {
        color0 = shortest_idx_v;
        if ( distance( checkerboard, color0 ) > 0.3) color0 = checkerboard;
            //color0 = texture(t2, vec3(v_st, panit));
        //color0 *= 1.0 - (mod(u_timer, period) / period);
    //} else if ((u_timer > 37000.) && (u_timer < 68000.)) {
    //	color0 = mixedv2;
        //if ( distance( checkerboard, color0 ) > 0.3) color0 = checkerboard;
            //color0 = texture(t2, vec3(v_st, panit));
        //color0 *= 1.0 - (mod(u_timer, period) / period);
    } else if ((u_timer > 52000.) && (u_timer < 82000.)) {
        vec4 panc = mixedv2;
        //panc *= 1.5 - (mod(u_timer, period) / period);
        if ( distance( plasma, mixedv1 ) < abs(sin(u_timer*0.0001))) color0 = mixedv1;
            else color0 = mixedv2;
        //color0 = (plasma*0.35+panc*0.65);
        //color0 *= 1.5 - (mod(u_timer, period) / period);
    } else if ((u_timer > 82000.) && (u_timer < 97000.)) {
        //color0 = plasma*vec4(1.,0.2,1.,1.);
        //vec4 panc = texture(t2, vec3(v_st, panit));
        //color0 = (plasma*0.65+panc*0.15+checkerboard*0.1+twirl*.1)*vec4(1.,0.5,1.,1.);
        //color0 += shortest_idx_v*0.5;
                        
        color0 = swirl * shortest_idx_v2 * vec4(1.9, 1.2, 1.7, 1.8);
        
    } else if ((u_timer > 97000.) && (u_timer < 126000.)) {
        color0 = rimina;
        //if ( distance( vec4(0.,0.,0.,1.0), color0 ) > 0.3) {
            color0 *= shortest_idx_v;
        //}
    } else if ((u_timer > 126000.) && (u_timer < 156000.)) {
        color0 = checkerboard;
        color0 *= 1.0 - (mod(u_timer, period) / period);
    } else if ((u_timer > 156000.) && (u_timer < 160100.)) {
        color0 = vec4(0.);
    }
*/
    
    color0 = vec4( c1 * color1 + c2 * color1 - twirl.xyz, 1. );

    if (u_timer < 10200.) { 
        // yellow concentric heart fade in
         color0 = vec4(c2 * color2 * easeIn(fadeIn)*0.25, 1.);
    } 
    else if ((u_timer > 10200.) && (u_timer < T_BRIDGE)) {
        // c1 pulsing red heart, c2 yellow hearts
        color0 = vec4( c1 * color1 * heartPulse*0.75/*pow(pow(sin(u_timer*0.00785)*0.5+0.5,2.),2.)*/ + c2 * color2, 1. );
    } 
    else if ((u_timer > T_BRIDGE) && (u_timer < T_BRIDGE_END/*31700*/)) {
        // tension building bridge, time slows, fade to black
        //color0 = vec4( /*c1 * color1 * heartPulse*/ 0./*pow(pow(sin(u_timer*0.00785)*0.5+0.5,2.),2.)*/ + c2 * color2, 1. );
        float tween = (u_timer - T_BRIDGE) / (T_BRIDGE_END - T_BRIDGE);
        color0 = vec4( c1 * color1 * (1.0-tween*tween*tween) + c2 * color2 * (1.0-tween) - tween*twirl.zzz, 1. );
        //color0 = vec4( /*1 * color2*.5 + c2 * color1*.5 -*/ twirl.zzz * sin(iTime), 1. );
    } 
    else if ((u_timer > T_BRIDGE_END) && (u_timer < 42090.)) {
        // heart beat red and red
        float rows = 512.0;
        float row = floor((v_st.y+iTime*0.01)*rows);// / rows;
        float interlace = mod(row, 2.0) > 0.5 ? 1.25 : 0.5;
        color0 = vec4( c1 * color2*interlace * (0.25+0.75*heartPulse) + c2 * (color1*0.5 + vec3(0.,0.,.05)) /* + twirl.xyz*0.5*/, 1. );
    }
    else if ((u_timer > 42090.) && (u_timer < 77450.)) {
        // bring in twirl
        float twirlIn = trange(42090., 46300.);
        if (u_timer > 50500.)
            twirlIn = 1.0-(0.25+0.75*trange(50500., 77450.));
        float tt = max(0., 1.0-(u_timer - 42090.)/(1000.*8.*0.860625));
        if (u_timer > 46300.)
            tt = max(tt, 1.0-(u_timer - 46300.)/(1000.*4.*0.860625));
        if (u_timer > 50500.)
            tt = max(tt, 1.0-(u_timer - 50500.)/(1000.*1.*0.860625));
        if (u_timer > 51300.)
            tt = max(tt, 1.0-(u_timer - 51300.)/(1000.*6.*0.860625));
        if (u_timer > 55816.)
            tt = max(tt, 1.0-(u_timer - 55816.)/(1000.*7.*0.860625));
        if (u_timer > 60945.)
            tt = max(tt, 1.0-(u_timer - 60945.)/(1000.*8.*0.860625));
        if (u_timer > 66000.)
            tt = max(tt, 1.0-(u_timer - 66000.)/(1000.*8.*0.860625));
        if (u_timer > 71300.)
            tt = max(tt, 1.0-(u_timer - 71300.)/(1000.*8.*0.860625));
        //color0 = vec4(tt);
        tt *= tt;
        //color0 = vec4( c1 * mix(color1+vec3(0.,pinch*0.125,pinch*0.5), color2, tt) * twirl.xyz + c2 * color1* (1.-heartPulse*.25) - twirl.xyz * heartPulse, 1. )*2.;
        color0 = vec4( c1 * (mix(vec3(0.), color2, tt)+vec3(0.,pinch*0.125,pinch*0.5))  + c2 * color1*0.75  - twirl.xyz *c2*0.75*twirlIn, 1. )*1.5;
        //color0 = twirl.xyzz;
        //color0.z += pinch;
        //+vec3(0.,0.,1.0)
    }
    else { // if (u_timer < 116000./*133720.*/) {
        //color0 = vec4( c1 * (color1+vec3(0.,pinch*0.125,pinch*0.5))  + c2 * color1*0.75  - twirl.xyz *c2*0.75, 1. )*2.;
        if (u_timer > 103310.)
            heartPulse = 1.;
        color0 = vec4( c1 * color2 * twirl.xyz + c2 * color1* (1.-heartPulse*.25) - twirl.xyz * heartPulse, 1. )*1.;
        if (u_timer > 104000.) {
            float tt = trange(104000., 130000.);
            //tt *= tt;
            color0 = mix(color0, vec4( c1 * mix(color2, color1, tt) * twirl.xyz + c2 * color1 * twirl.xyz + twirl.xxx*max(0.,1.-c1-c2)*0.5, 1. )*2., tt);
        }
    }

    
    if (u_timer > 120000.) { // fade to black at end
        color0 *= 1.0-fadeOut;
    }

    color = color0;
}


