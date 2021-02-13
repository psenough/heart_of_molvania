#version 300 es

precision highp float;
precision highp int;

uniform int layer;
uniform float u_timer;
uniform float u_delta_timer;
uniform vec2 u_resolution;

in vec2 v_st;

out vec4 color;
float iTime;
vec2 uv;
float PI = 3.14159265359;

vec2 RotCS(vec2 p, float c, float s) {
    return vec2( p.x*c+p.y*s,
                -p.x*s+p.y*c);
}

//iq 2d simplex noise
vec2 hash(vec2 p) {
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

float HASHSCALE1 = 443.8975;
float hash13(vec3 p3){
    p3  = fract(p3 * HASHSCALE1);
    p3 += dot(p3, p3.yzx + 19.19);
    return fract((p3.x + p3.y) * p3.z);
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
float trange(float a, float b) {
    return max(0., min(1., (u_timer - a) / (b - a)));
}

void main() {
    vec4 color0 = vec4(0.);
    
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
        float rows = mix(32., 256., tween);
        float amp = 0.125 * ramp;
        float row = floor(v_st.y*rows) / rows;
        uv.x += hash13(vec3(row, 0., 0.)) * amp;
    }
    else if (u_timer > 77600.) { // vert row distort
        float rows = 256.;
        float amp = max(0., sin(u_timer*0.00785))*0.075;
        float row = floor(v_st.x*rows) / rows;
        vec2 hh = hash(vec2(u_timer/1000000., row));
        if (hh.x > 0.9)
            uv.y += hh.y * amp;
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
        iTime += mix(0., 4096.*heartZoomRate, mx);
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

    color0 = vec4( c1 * color1 + c2 * color1 - twirl.xyz, 1. );

    if (u_timer < 10200.) { 
        // yellow concentric heart fade in
         color0 = vec4(c2 * color2 * easeIn(fadeIn)*0.25, 1.);
    } 
    else if ((u_timer > 10200.) && (u_timer < T_BRIDGE)) {
        // c1 pulsing red heart, c2 yellow hearts
        color0 = vec4( c1 * color1 * heartPulse*0.75 + c2 * color2, 1. );
    } 
    else if ((u_timer > T_BRIDGE) && (u_timer < T_BRIDGE_END/*31700*/)) {
        // tension building bridge, time slows, fade to black
        float tween = (u_timer - T_BRIDGE) / (T_BRIDGE_END - T_BRIDGE);
        color0 = vec4( c1 * color1 * (1.0-tween*tween*tween) + c2 * color2 * (1.0-tween) - tween*twirl.zzz, 1. );
    } 
    else if ((u_timer > T_BRIDGE_END) && (u_timer < 42090.)) {
        // heart beat red and red
        float rows = 512.0;
        float row = floor((v_st.y+iTime*0.01)*rows);// / rows;
        float interlace = mod(row, 2.0) > 0.5 ? 1.25 : 0.5;
        color0 = vec4( c1 * color2*interlace * (0.25+0.75*heartPulse) + c2 * (color1*0.5 + vec3(0.,0.,.05)), 1. );
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
        tt *= tt;
        color0 = vec4( c1 * (mix(vec3(0.), color2, tt)+vec3(0.,pinch*0.125,pinch*0.5))  + c2 * color1*0.75  - twirl.xyz *c2*0.75*twirlIn, 1. )*1.5;
    }
    else {
        if (u_timer > 103310.)
            heartPulse = 1.;
        color0 = vec4( c1 * color2 * twirl.xyz + c2 * color1* (1.-heartPulse*.25) - twirl.xyz * heartPulse, 1. )*1.;
        if (u_timer > 104000.) {
            float tt = trange(104000., 130000.);
            color0 = mix(color0, vec4( c1 * mix(color2, color1, tt) * twirl.xyz + c2 * color1 * twirl.xyz + twirl.xxx*max(0.,1.-c1-c2)*0.5, 1. )*2., tt);
        }
    }
    
    if (u_timer > 120000.) { // fade to black at end
        color0 *= 1.0-fadeOut;
    }

    color = color0;
}
