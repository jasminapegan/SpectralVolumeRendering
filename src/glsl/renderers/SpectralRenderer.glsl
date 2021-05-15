// #package glsl/shaders

// #include ../mixins/SpectralPhoton.glsl
// #include ../mixins/rand.glsl
// #include ../mixins/unprojectRand.glsl
// #include ../mixins/intersectCube.glsl
// #include ../mixins/spectrumHelpers.glsl

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

uniform float uPhotonFreq;
uniform mat4 uSpectrum;

in vec2 vPosition;

layout (location = 0) out vec4 oPosition;
layout (location = 1) out vec4 oDirection;
layout (location = 2) out vec4 oTransmittance;
layout (location = 3) out vec4 oRadiance;

@rand
@unprojectRand
@intersectCube
@spectrumHelpers

void resetPhoton(inout vec2 randState, inout SpectralPhoton photon) {
    vec3 from, to;
    unprojectRand(randState, vPosition, uMvpInverseMatrix, uInverseResolution, uBlur, from, to);
    photon.direction = normalize(to - from);
    photon.bounces = 0u;
    vec2 tbounds = max(intersectCube(from, photon.direction), 0.0);
    photon.position = from + tbounds.x * photon.direction;
    photon.transmittance = vec3(1);
    photon.frequency = uPhotonFreq; //sampleFrequency(randState);  // get random from spectrum
    photon.spectrum = float[8](uSpectrum[0][0], uSpectrum[0][1], uSpectrum[0][2], uSpectrum[0][3],
                                uSpectrum[1][0], uSpectrum[1][1], uSpectrum[1][2], uSpectrum[1][3]);
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
    photon.spectrum = float[8](uSpectrum[0][0], uSpectrum[0][1], uSpectrum[0][2], uSpectrum[0][3],
                                uSpectrum[1][0], uSpectrum[1][1], uSpectrum[1][2], uSpectrum[1][3]);

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
        uint j = idx;

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
            //photon.spectrum[idx] *= 1.0 - length(dx) * volumeSample.x;  // assuming only this band is absorbed

        } else if (r.y < PAbsorption) {
            // absorption
            float weightA = muAbsorption / (uMajorant * PAbsorption);
            photon.transmittance *= 1.0 - weightA;
            // spectrum -- update
            //photon.spectrum[idx] *= 1.0 - length(dx) * volumeSample.x;  // assuming only this band is absorbed

        } else if (r.y < PAbsorption + PScattering) {
            // scattering
            r = rand(r);
            float weightS = muScattering / (uMajorant * PScattering);
            photon.transmittance *= volumeSample.rgb * weightS;
            photon.direction = sampleHenyeyGreenstein(uScatteringBias, r, photon.direction);
            photon.bounces++;
            // spectrum -- update
            float s = photon.spectrum[idx] * 1.0 - float(length(dx)) * volumeSample.x;
            //photon.spectrum[idx] = s; //*= 1.0 - float(length(dx)) * volumeSample.x;  // assuming only this band is absorbed
            photon.spectrum[j] = 0.0;

        } else {
            // null collision
            float weightN = muNull / (uMajorant * PNull);
            photon.transmittance *= weightN;
        }
    }
    oPosition = vec4(photon.position, 0);
    oDirection = vec4(photon.direction, float(photon.bounces));
    //oTransmittance = vec4(photon.transmittance, 0); // color actually
    oTransmittance = vec4(spectrumToRgb(photon.spectrum), 0.0);
    /*vec4(photon.spectrum[0], photon.spectrum[1], photon.spectrum[2], photon.spectrum[3]), 
    vec4(photon.spectrum[4], photon.spectrum[5], photon.spectrum[6], photon.spectrum[7])),
    0); //0.0, 0.0, 0.0, 0.0);*/
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

#define EPS 1e-5

@SpectralPhoton

uniform mat4 uMvpInverseMatrix;
uniform vec2 uInverseResolution;
uniform float uRandSeed;
uniform float uBlur;
uniform float uPhotonFreq;
uniform mat4 uSpectrum;

in vec2 vPosition;

layout (location = 0) out vec4 oPosition;
layout (location = 1) out vec4 oDirection;
layout (location = 2) out vec4 oTransmittance;
layout (location = 3) out vec4 oRadiance;

@rand
@unprojectRand
@intersectCube
@spectrumHelpers

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
    photon.frequency = uPhotonFreq;
    photon.bounces = 0u;
    photon.samples = 0u;
    photon.spectrum = float[8](uSpectrum[0][0], uSpectrum[0][1], uSpectrum[0][2], uSpectrum[0][3],
                                uSpectrum[1][0], uSpectrum[1][1], uSpectrum[1][2], uSpectrum[1][3]);
    oPosition = vec4(photon.position, 0);
    oDirection = vec4(photon.direction, float(photon.bounces));
    //oTransmittance = vec4(photon.transmittance, 0);
    oTransmittance = vec4(spectrumToRgb(photon.spectrum), 0.0);
    /*vec4(photon.spectrum[0], photon.spectrum[1], photon.spectrum[2], photon.spectrum[3]), 
    vec4(photon.spectrum[4], photon.spectrum[5], photon.spectrum[6], photon.spectrum[7]),
    0));*/
    oRadiance = vec4(photon.radiance, float(photon.samples));
}
