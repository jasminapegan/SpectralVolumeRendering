// #package glsl/mixins

// #section spectrumHelpers


const uint len = 16u;

vec3 frequencyToRGB(float freq) {
    float Wavelength = pow(10.0, 9.0) * 299792458.0 / freq;
    float Gamma = 0.80;
    float IntensityMax = 255.0;
    float factor;
    float Red, Green, Blue;

    if((Wavelength >= 380.0) && (Wavelength < 440.0)) {
        Red = -(Wavelength - 440.0) / (440.0 - 380.0);
        Green = 0.0;
        Blue = 1.0;
    } else if((Wavelength >= 440.0) && (Wavelength < 490.0)) {
        Red = 0.0;
        Green = (Wavelength - 440.0) / (490.0 - 440.0);
        Blue = 1.0;
    } else if((Wavelength >= 490.0) && (Wavelength < 510.0)) {
        Red = 0.0;
        Green = 1.0;
        Blue = -(Wavelength - 510.0) / (510.0 - 490.0);
    } else if((Wavelength >= 510.0) && (Wavelength < 580.0)) {
        Red = (Wavelength - 510.0) / (580.0 - 510.0);
        Green = 1.0;
        Blue = 0.0;
    } else if((Wavelength >= 580.0) && (Wavelength < 645.0)) {
        Red = 1.0;
        Green = -(Wavelength - 645.0) / (645.0 - 580.0);
        Blue = 0.0;
    } else if((Wavelength >= 645.0) && (Wavelength < 781.0)) {
        Red = 1.0;
        Green = 0.0;
        Blue = 0.0;
    } else {
        Red = 0.0;
        Green = 0.0;
        Blue = 0.0;
    }

    // Let the intensity fall off near the vision limits; this is line 155
    if((Wavelength >= 380.0) && (Wavelength < 420.0)) {
        factor = 0.3 + 0.7 * (Wavelength - 380.0) / (420.0 - 380.0);
    } else if((Wavelength >= 420.0) && (Wavelength < 701.0)) {
        factor = 1.0;
    } else if((Wavelength >= 701.0) && (Wavelength < 781.0)) {
        factor = 0.3 + 0.7 * (780.0 - Wavelength) / (780.0 - 700.0);
    } else {
        factor = 0.0;
    }

    // Don't want 0^x = 1 for x <> 0
    int r = Red == 0.0 ? 0 : int(round(IntensityMax * pow(Red * factor, Gamma)));
    int g = Green == 0.0 ? 0 : int(round(IntensityMax * pow(Green * factor, Gamma)));
    int b = Blue == 0.0 ? 0 : int(round(IntensityMax * pow(Blue * factor, Gamma)));

    return vec3(r, g, b);
}

vec3 spectrumToRgb(float spectrum[len]) {
    vec3 rgb = vec3(0.0, 0.0, 0.0);
    for (uint i=0u; i<len; i++) {
        rgb += frequencyToRGB(spectrum[i]);
    }
    return rgb;
}

vec3 rgbToHsv(vec3 c)
{
    vec4 K = vec4(0.0, -1.0 / 3.0, 2.0 / 3.0, -1.0);
    vec4 p = mix(vec4(c.bg, K.wz), vec4(c.gb, K.xy), step(c.b, c.g));
    vec4 q = mix(vec4(p.xyw, c.r), vec4(c.r, p.yzx), step(p.x, c.r));

    float d = q.x - min(q.w, q.y);
    return vec3(abs(q.z + (q.w - q.y) / (6.0 * d + EPS)), d / (q.x + EPS), q.x);
}

uint rgbToSpectrumIndex(vec3 rgb) {
    float hue = rgbToHsv(rgb)[0];
    float wavelength = 620.0 - 170.0 / 270.0 * hue;
    float freq = 299792458.0 / (wavelength + EPS);
    // freq = 4 + 4 * i / len;
    // (freq - 4) * len = 4 * i
    // i = (freq - 4) * len/4
    return uint((freq - 4.0) * float(len) / 4.0);
}


vec3 limitFreq(float freq, vec3 color) {
    vec3 freqColor = frequencyToRGB(freq);
    return vec3(min(freqColor.x, color.x), min(freqColor.y, color.y), min(freqColor.z, color.z));
}


vec3 cropColor(vec3 color, vec3 filterColor, float coef) {
    vec3 ratio = vec3(
        max(0.0, color.x - filterColor.x) / color.x,
        max(0.0, color.y - filterColor.y) / color.y,
        max(0.0, color.z - filterColor.z) / color.z
    );
    return vec3(
        color.x * coef * ratio.x,
        color.y * coef * ratio.y,
        color.z * coef * ratio.z
    );
}
