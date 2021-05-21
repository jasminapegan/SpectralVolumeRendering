// #package js/main

// #include ../AbstractDialog.js
// #include ../../TransferFunctionWidget.js
// #include ../../AbsorptionWidget.js

// #include ../../../uispecs/renderers/SpectralRendererDialog.json

class SpectralRendererDialog extends AbstractDialog {

constructor(renderer, options) {
    super(UISPECS.SpectralRendererDialog, options);

    this._renderer = renderer;

    this._handleChange = this._handleChange.bind(this);
    this._handleTFChange = this._handleTFChange.bind(this);
    this._handleAChange = this._handleAChange.bind(this);

    this._binds.extinction.addEventListener('input', this._handleChange);
    this._binds.albedo.addEventListener('change', this._handleChange);
    this._binds.bias.addEventListener('change', this._handleChange);
    this._binds.ratio.addEventListener('change', this._handleChange);
    this._binds.bounces.addEventListener('input', this._handleChange);
    this._binds.steps.addEventListener('input', this._handleChange);
    this._binds.n.addEventListener('input', this._handleChange);
    this._binds.temperature.addEventListener('input', this._handleChange);

    this._tfwidget = new TransferFunctionWidget();
    this._binds.tfcontainer.add(this._tfwidget);
    this._tfwidget.addEventListener('change', this._handleTFChange);

    this._awidget = new AbsorptionWidget();
    this._binds.acontainer.add(this._awidget);
    this._awidget.addEventListener('change', this._handleAChange);
}

destroy() {
    this._tfwidget.destroy();
    this._awidget.destroy();
    super.destroy();
}

_handleChange() {
    const extinction = this._binds.extinction.getValue();
    const albedo     = this._binds.albedo.getValue();
    const bias       = this._binds.bias.getValue();
    const ratio      = this._binds.ratio.getValue();
    const bounces    = this._binds.bounces.getValue();
    const steps      = this._binds.steps.getValue();
    const n          = this._binds.n.getValue();
    const temperature = this._binds.temperature.getValue();

    this._renderer.absorptionCoefficient = extinction * (1 - albedo);
    this._renderer.scatteringCoefficient = extinction * albedo;
    this._renderer.scatteringBias = bias;
    this._renderer.majorant = extinction * ratio;
    this._renderer.maxBounces = bounces;
    this._renderer.steps = steps;
    this._renderer.n = n;
    this._renderer.temperature = temperature;

    this._renderer.reset();
}

_handleTFChange() {
    this._renderer.setTransferFunction(this._tfwidget.getTransferFunction());
    this._renderer.reset();
}

_handleAChange() {
    this._renderer.setAbsorption(this._awidget.getAbsorption());
    this._renderer.reset();
}

}
