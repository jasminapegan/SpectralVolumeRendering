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
uniform uint uIdx;

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
    
    float Wavelength = 299792458.0 / photon.frequency;
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

    photon.transmittance = vec3(1); //vec3(r, g, b);
    photon.frequency = uPhotonFreq;
}

vec4 sampleEnvironmentMap(vec3 d) {
    vec2 texCoord = vec2(atan(d.x, -d.z), asin(-d.y) * 2.0) * M_INVPI * 0.5 + 0.5;
    return texture(uEnvironment, texCoord);
}

vec4 sampleVolumeColor(vec3 position) {
    vec2 volumeSample = texture(uVolume, position).rg;
    vec4 transferSample = texture(uTransferFunction, volumeSample);
    return transferSample;
}

vec4 sampleAbsorptionColor(vec3 position) {  
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
    photon.frequency = uPhotonFreq;
    vec2 mappedPosition = vPosition * 0.5 + 0.5;
    photon.position = texture(uPosition, mappedPosition).xyz;
    vec4 directionAndBounces = texture(uDirection, mappedPosition);
    photon.direction = directionAndBounces.xyz;
    photon.bounces = uint(directionAndBounces.w + 0.5);
    photon.transmittance = texture(uTransmittance, mappedPosition).rgb;
    //photon.transmittance = limitFreq(photon.frequency, photon.transmittance);
    vec4 radianceAndSamples = texture(uRadiance, mappedPosition);
    photon.radiance = radianceAndSamples.rgb;
    photon.samples = uint(radianceAndSamples.w + 0.5);

    vec2 r = rand(vPosition * uRandSeed);
    for (uint i = 0u; i < uSteps; i++) {
        r = rand(r);
        float t = -log(r.x) / uMajorant;
        vec3 dx = t * photon.direction;
        photon.position += dx;

        vec4 volumeSample = sampleVolumeColor(photon.position);
        float density = texture(uVolume, photon.position).r;
        vec4 absorptionSample = sampleAbsorptionColor(photon.position);
        float muAbsorption = density * absorptionSample.w;
        float muScattering = density * uScatteringCoefficient; //volumeSample.a * uScatteringCoefficient;
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
            photon.samples++;
            photon.radiance += (radiance - photon.radiance) / float(photon.samples);
            resetPhoton(r, photon);
        } else if (photon.bounces >= uMaxBounces) {
            // max bounces achieved -> only estimate transmittance
            float weightAS = (muAbsorption + muScattering) / uMajorant;
            photon.transmittance *= 1.0 - weightAS;
        }
        else if (idx == uIdx) {//r.y < PAbsorption) {
            // absorption
            //float weightA = muAbsorption / (uMajorant * PAbsorption);
            //photon.transmittance *= 1.0 - weightA;
            // if photon is in the same spectral band as absorbed color, execute absorption
            float factor = 1.0 - float(length(dx)) * muAbsorption;
            photon.transmittance *= volumeSample.rgb * factor;
        }
        else if (r.y < PAbsorption + PScattering) {
            // scattering
            if (uIdx != idx) {
                r = rand(r);
                float weightS = muScattering / (uMajorant * PScattering);
                //photon.transmittance *= volumeSample.rgb * weightS;
                vec3 factor = vec3(1) - float(length(dx)) * muAbsorption * absorptionSample.rgb;
                photon.transmittance *= volumeSample.rgb * factor;
                /*if (uIdx != idx) {
                    photon.transmittance *= volumeSample.a * absorptionSample.w;
                } else {
                    photon.transmittance *= factor;
                    //photon.transmittance = cropColor(photon.transmittance, absorptionSample.rgb, factor);
                }*/
                photon.direction = sampleHenyeyGreenstein(uScatteringBias, r, photon.direction);
                photon.bounces++;
            }
        }
        else {
            // null collision
            float weightN = muNull / (uMajorant * PNull);
            photon.transmittance *= weightN;
        }
    }

    oPosition = vec4(photon.position, 0);
    oDirection = vec4(photon.direction, float(photon.bounces));
    oTransmittance = vec4(photon.transmittance, 0);
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
uniform float uPhotonFreq;

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
    photon.frequency = uPhotonFreq;
    photon.bounces = 0u;
    photon.samples = 0u;
    oPosition = vec4(photon.position, 0);
    oDirection = vec4(photon.direction, float(photon.bounces));
    oTransmittance = vec4(photon.transmittance, 0);
    oRadiance = vec4(photon.radiance, float(photon.samples));
}
