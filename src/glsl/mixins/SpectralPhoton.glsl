// #package glsl/mixins

// #section SpectralPhoton

struct SpectralPhoton {
    vec3 position;
    vec3 direction;
    vec3 transmittance;
    vec3 radiance;
    uint bounces;
    uint samples;
    float frequency;  // 4.3-7.5
};
