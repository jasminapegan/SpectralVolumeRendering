/*// #package glsl/shaders

// #include ../mixins/SpectralPhoton.glsl
// #include ../mixins/rand.glsl
// #include ../mixins/unprojectRand.glsl
// #include ../mixins/intersectCube.glsl

// #section SpectralGenerate/vertex

void main() {}

// #section SpectralGenerate/fragment

void main() {}

// #section SpectralIntegrate/vertex

#version 300 es

layout (location = 0) in vec2 aPosition;

out vec2 vPosition;

void main() {
    vPosition = aPosition;
    gl_Position = vec4(aPosition, 0.0, 1.0);
}

// #section SpectralIntegrate/fragment

#version 300 es
precision mediump float;

#define M_INVPI 0.31830988618
#define M_2PI 6.28318530718
#define EPS 1e-5

@SpectralPhoton

uniform mediump sampler2D uPosition;
uniform mediump sampler2D uDirection;
uniform mediump sampler2D uTransmittance;
uniform mediump sampler2D uRadiance;

uniform mediump sampler3D uVolume;
uniform mediump sampler2D uTransferFunction;
uniform mediump sampler2D uAbsorption;
uniform mediump sampler2D uEnvironment;

uniform mat4 uMvpInverseMatrix;
uniform vec2 uInverseResolution;
uniform float uRandSeed;
uniform float uBlur;

uniform float uAbsorptionCoefficient;
uniform float uScatteringCoefficient;
uniform float uScatteringBias;
uniform float uMajorant;
uniform uint uMaxBounces;
uniform uint uSteps;

in vec2 vPosition;

layout (location = 0) out vec4 oPosition;
layout (location = 1) out vec4 oDirection;
layout (location = 2) out vec4 oTransmittance;
layout (location = 3) out vec4 oRadiance;

@rand
@unprojectRand
@intersectCube

float sampleFrequency(vec2 randState) {
    float[] spectrum = float[8](2.7, 3.3, 2.6, 2.5, 2.3, 2.1, 2.0, 1.75);
    const int len = 8;
    float sum = 0.0;
    for (int i=0; i < len; i++) {
        sum += spectrum[i];
    }
    float[] cumulativeSpectrum = float[len](0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    cumulativeSpectrum[0] = spectrum[0] / sum;
    for (int i=1; i < len; i++) {
        cumulativeSpectrum[i] = (cumulativeSpectrum[i-1] + spectrum[i]) / sum;
    }

    float u = rand(randState)[0]; //random iz [0,1]? preveri ce je res
    int i = 0;
    while (cumulativeSpectrum[i] > u) {
        i++;
    }
    return 4.0 + float(4 * i) / float(len); // return actual freq between 4 and 7.5
}

//Taken from Earl F. Glynn's web page:
// <a href="http://www.efg2.com/Lab/ScienceAndEngineering/Spectra.htm">Spectra Lab Report</a>
vec3 frequencyToRGB(float freq) {
    float Wavelength = 299792458.0 / freq;
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

vec3 spectrumToRgb(float[8] spectrum) {
    vec3 rgb = vec3(0.0, 0.0, 0.0);
    for (uint i=0u; i<8u; i++) {
        rgb += frequencyToRGB(spectrum[i]);
    }
    return rgb;
}

void resetPhoton(inout vec2 randState, inout SpectralPhoton photon) {
    vec3 from, to;
    unprojectRand(randState, vPosition, uMvpInverseMatrix, uInverseResolution, uBlur, from, to);
    photon.direction = normalize(to - from);
    photon.bounces = 0u;
    vec2 tbounds = max(intersectCube(from, photon.direction), 0.0);
    photon.position = from + tbounds.x * photon.direction;
    photon.transmittance = vec3(1);
    photon.frequency = sampleFrequency(randState);  // get random from spectrum
}

vec4 sampleEnvironmentMap(vec3 d) {
    vec2 texCoord = vec2(atan(d.x, -d.z), asin(-d.y) * 2.0) * M_INVPI * 0.5 + 0.5;
    return texture(uEnvironment, texCoord);
}

vec4 sampleVolumeColor(vec3 position) {    // nekje tle dodaj kolicnik za lom
    vec2 volumeSample = texture(uVolume, position).rg;
    vec4 transferSample = texture(uTransferFunction, volumeSample);
    return transferSample;
}

vec4 sampleAbsorptionColor(vec3 position) {    // nekje tle dodaj kolicnik za lom
    vec2 volumeSample = texture(uVolume, position).rg;
    vec4 absorptionSample = texture(uAbsorption, volumeSample);
    return absorptionSample;
}

vec3 randomDirection(vec2 U) {
    float phi = U.x * M_2PI;
    float z = U.y * 2.0 - 1.0;
    float k = sqrt(1.0 - z * z);
    return vec3(k * cos(phi), k * sin(phi), z);
}

float sampleHenyeyGreensteinAngleCosine(float g, float U) {
    float g2 = g * g;
    float c = (1.0 - g2) / (1.0 - g + 2.0 * g * U);
    return (1.0 + g2 - c * c) / (2.0 * g);
}

vec3 sampleHenyeyGreenstein(float g, vec2 U, vec3 direction) {
    // generate random direction and adjust it so that the angle is HG-sampled
    vec3 u = randomDirection(U);
    if (abs(g) < EPS) {
        return u;
    }
    float hgcos = sampleHenyeyGreensteinAngleCosine(g, fract(sin(U.x * 12345.6789) + 0.816723));
    float lambda = hgcos - dot(direction, u);
    return normalize(u + lambda * direction);
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
    float[] spectrum = float[8](2.7, 3.3, 2.6, 2.5, 2.3, 2.1, 2.0, 1.75);
    float hue = rgbToHsv(rgb)[0];
    float wavelength = 620.0 - 170.0 / 270.0 * hue;
    float freq = 299792458.0 / wavelength;
    uint i = 0u;
    while (spectrum[i] < freq) {
        i++;
    }
    return i;
}

void main() {
    SpectralPhoton photon;
    vec2 mappedPosition = vPosition * 0.5 + 0.5;
    photon.position = texture(uPosition, mappedPosition).xyz;
    vec4 directionAndBounces = texture(uDirection, mappedPosition);
    photon.direction = directionAndBounces.xyz;
    photon.bounces = uint(directionAndBounces.w + 0.5);
    photon.transmittance = texture(uTransmittance, mappedPosition).rgb;
    vec4 radianceAndSamples = texture(uRadiance, mappedPosition);
    photon.radiance = radianceAndSamples.rgb;
    photon.samples = uint(radianceAndSamples.w + 0.5);
    photon.spectrum = float[8](0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

    vec2 r = rand(vPosition * uRandSeed);
    for (uint i = 0u; i < uSteps; i++) {
        r = rand(r);
        float t = -log(r.x) / uMajorant;
        vec3 dx = t * photon.direction;
        photon.position += dx;

        vec4 volumeSample = sampleVolumeColor(photon.position);
        vec4 absorptionSample = sampleAbsorptionColor(photon.position);
        float muAbsorption = volumeSample.a * uAbsorptionCoefficient;
        float muScattering = volumeSample.a * uScatteringCoefficient;
        float muNull = uMajorant - muAbsorption - muScattering;
        float muMajorant = muAbsorption + muScattering + abs(muNull);
        float PNull = abs(muNull) / muMajorant;
        float PAbsorption = muAbsorption / muMajorant;
        float PScattering = muScattering / muMajorant;
        uint idx = rgbToSpectrumIndex(absorptionSample.rgb);  // index of absorbed band

        if (any(greaterThan(photon.position, vec3(1))) || any(lessThan(photon.position, vec3(0)))) {
            // out of bounds
            vec4 envSample = sampleEnvironmentMap(photon.direction);
            vec3 radiance = photon.transmittance * envSample.rgb;
            photon.radiance += (radiance - photon.radiance) / float(photon.samples);
            photon.samples++;
            resetPhoton(r, photon);

        } else if (photon.bounces >= uMaxBounces) {
            // max bounces achieved -> only estimate transmittance
            float weightAS = (muAbsorption + muScattering) / uMajorant;
            photon.transmittance *= 1.0 - weightAS;
            // spectrum -- update
            photon.spectrum[idx] *= 1.0 - length(dx) * volumeSample.a;  // assuming only this band is absorbed

        } else if (r.y < PAbsorption) {
            // absorption
            float weightA = muAbsorption / (uMajorant * PAbsorption);
            photon.transmittance *= 1.0 - weightA;
            // spectrum -- update
            photon.spectrum[idx] *= 1.0 - length(dx) * volumeSample.a;  // assuming only this band is absorbed

        } else if (r.y < PAbsorption + PScattering) {
            // scattering
            r = rand(r);
            float weightS = muScattering / (uMajorant * PScattering);
            photon.transmittance *= volumeSample.rgb * weightS;
            photon.direction = sampleHenyeyGreenstein(uScatteringBias, r, photon.direction);
            photon.bounces++;
            // spectrum -- update
            photon.spectrum[idx] *= 1.0 - length(dx) * volumeSample.a;  // assuming only this band is absorbed

        } else {
            // null collision
            float weightN = muNull / (uMajorant * PNull);
            photon.transmittance *= weightN;
        }
    }
    oPosition = vec4(photon.position, 0);
    oDirection = vec4(photon.direction, float(photon.bounces));
    //oTransmittance = vec4(photon.transmittance, 0); // color actually
    oTransmittance = vec4(spectrumToRgb(photon.spectrum), 0);
    oRadiance = vec4(photon.radiance, float(photon.samples));
}

// #section SpectralRender/vertex

#version 300 es

layout (location = 0) in vec2 aPosition;
out vec2 vPosition;

void main() {
    vPosition = (aPosition + 1.0) * 0.5;
    gl_Position = vec4(aPosition, 0.0, 1.0);
}

// #section SpectralRender/fragment

#version 300 es
precision mediump float;

uniform mediump sampler2D uColor;

in vec2 vPosition;
out vec4 oColor;

void main() {
    oColor = vec4(texture(uColor, vPosition).rgb, 1);
}

// #section SpectralReset/vertex

#version 300 es

layout (location = 0) in vec2 aPosition;

out vec2 vPosition;

void main() {
    vPosition = aPosition;
    gl_Position = vec4(aPosition, 0.0, 1.0);
}

// #section SpectralReset/fragment

#version 300 es
precision mediump float;

@SpectralPhoton

uniform mat4 uMvpInverseMatrix;
uniform vec2 uInverseResolution;
uniform float uRandSeed;
uniform float uBlur;

in vec2 vPosition;

layout (location = 0) out vec4 oPosition;
layout (location = 1) out vec4 oDirection;
layout (location = 2) out vec4 oTransmittance;
layout (location = 3) out vec4 oRadiance;

@rand
@unprojectRand
@intersectCube

void main() {
    SpectralPhoton photon;
    vec3 from, to;
    vec2 randState = rand(vPosition * uRandSeed);
    unprojectRand(randState, vPosition, uMvpInverseMatrix, uInverseResolution, uBlur, from, to);
    photon.direction = normalize(to - from);
    vec2 tbounds = max(intersectCube(from, photon.direction), 0.0);
    photon.position = from + tbounds.x * photon.direction;
    photon.transmittance = vec3(1);
    photon.radiance = vec3(1);
    photon.bounces = 0u;
    photon.samples = 0u;
    photon.spectrum = photon.spectrum;
    oPosition = vec4(photon.position, 0);
    oDirection = vec4(photon.direction, float(photon.bounces));
    //oTransmittance = vec4(photon.transmittance, 0);
    oTransmittance = vec4(spectrumToRgb(photon.spectrum), 0);
    oRadiance = vec4(photon.radiance, float(photon.samples));
}
*/