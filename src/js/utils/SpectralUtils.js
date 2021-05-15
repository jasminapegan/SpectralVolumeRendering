// #package js/main

class SpectralUtils {


static sampleFrequency() {
    var spectrumWeights = [2.7, 3.3, 2.6, 2.5, 2.3, 2.1, 2.0, 1.75];
    var len = 8;

    var sum = 0.0;
    for (var i=0; i < len; i++) {
        sum += spectrumWeights[i];
    }
    var cumulativeSpectrum = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    
    for (var i=1; i < len; i++) {
        cumulativeSpectrum[i] = (cumulativeSpectrum[i-1] + spectrumWeights[i]) / sum;
    }

    var u = Math.random(); //random iz [0,1]? preveri ce je res
    var i = 0;
    while (cumulativeSpectrum[i] > u) {
        i++;
    }
    return 4 + 4 * i / len; // return actual freq between 4 and 7.5
}

}