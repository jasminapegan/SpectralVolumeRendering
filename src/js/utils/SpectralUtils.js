// #package js/main

class SpectralUtils {

static blackBodyWeight(lambda, dist, T) {
    // dist is freq category size
    var c1 = 3.74183 * Math.pow(10, -16); //2 pi h c^2
    var c2 = 0.014388;  // h c / kB
    var e = 2.71828;
    return c1 * Math.pow(lambda, -5) * dist / (Math.pow(e, c2 / (lambda * T)) - 1);
}

static blackBodyWeights(n, minLambda, maxLambda, T) {
    var weights = new Array(n);
    var dist = (maxLambda - minLambda) / n;
    for (var i=0; i < n; i++) {
        weights[i] = this.blackBodyWeight(minLambda + i * dist + dist/2, dist, T);
    }
    return weights;
}

static cumulativeSpectrum(n, minLambda, maxLambda, T)
{
    var spectrum = this.blackBodyWeights(n, minLambda, maxLambda, T) ;
    var num = 16;
    var minLambda = 400 * Math.pow(10, -9);
    var maxLambda = 700 * Math.pow(10, -9);
    var spectrumWeights = this.blackBodyWeights(num, minLambda, maxLambda, T);

    var sum = 0.0;
    for (var i=0; i < num; i++) {
        sum += spectrumWeights[i];
    }
    var cumulativeSpectrum = [0.0];
    for (var i=1; i < num; i++) {
        cumulativeSpectrum[i] = (cumulativeSpectrum[i-1] + spectrumWeights[i]);
    }
    return [cumulativeSpectrum, sum];
}

static sampleFrequency(T, spec, sum) {
    var c = 299792458.0;
    var num = 16;
    var u = Math.random() * sum;
    var i = 0;
    var minLambda = 400 * Math.pow(10, -9);
    var maxLambda = 700 * Math.pow(10, -9);
    while (spec[i] < u) {
        i++;
    }
    var lambda = minLambda + (maxLambda - minLambda) * i / num;
    return [c / lambda, i];
}

}