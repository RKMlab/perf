/*-- PIE PLOT -------------------------------------------------------------------------------*/
const pieDatum = function(data, type, repeats, group) {
    let datum;
    let pieData = {};
    for (let d in data) {
        if (repeats.indexOf(d) > -1) {
            pieData[d] = data[d];
        }
    }
    let totalFreq = 0;
    let totalBases = 0;
    for (let i in pieData) {
        totalFreq += _.sum(pieData[i]);
        totalBases += (_.sum(pieData[i])) * (i.length);
    }

    if (group == 1) {
        let kmerDataObj = {};
        for (let k in kMers) { kmerDataObj[kMers[k]] = 0; }
        datum = _.map(Object.keys(pieData), function(d) {
            if (type == 'Frequency') { kmerDataObj[kmerObj[d.length]] += _.sum(pieData[d]); } else { kmerDataObj[kmerObj[d.length]] += _.sum(pieData[d]) * d.length; }
        });
        datum = [];
        for (let k in kmerDataObj) { datum.push({ name: k, y: kmerDataObj[k] }); }
    } else {
        datum = _.map(Object.keys(pieData), function(d) {
            let point;
            if (type == 'Frequency') { point = { name: d, y: _.sum(pieData[d]) } } else { point = { name: d, y: _.sum(pieData[d]) * d.length }; }
            return point;
        });
    }
    //  console.log(datum)
    return datum
}

const piePlot = function(repeats) {

    const dataTypeRadios = document.getElementsByName('data-type');
    console.log(dataTypeRadios)
    let dataType;
    for (let i in dataTypeRadios) {
        let radio = dataTypeRadios[i];
        if (radio.tagName == 'INPUT') {
            if (radio.checked) {
                dataType = radio.value;
            }
        }
    }
    console.log(dataType);
    const lenGroupRadios = document.getElementsByName('len-group');
    let lenGroup;
    for (let i in lenGroupRadios) {
        let radio = lenGroupRadios[i];
        if (radio.tagName == 'INPUT') {
            if (radio.checked) {
                lenGroup = radio.value;
            }
        }
    }
    let plotData = pieDatum(plotInfo['len'], dataType, repeats, lenGroup);
    Highcharts.chart('pie-plot-svg', {
        chart: {
            plotBackgroundColor: null,
            plotBorderWidth: null,
            plotShadow: false,
            type: 'pie',
            marginTop: 40
        },
        title: {
            text: null
        },
        tooltip: {
            pointFormat: '{series.name}: <b>{point.percentage:.1f}%</b>'
        },
        plotOptions: {
            pie: {
                allowPointSelect: true,
                cursor: 'pointer',
                dataLabels: {
                    enabled: false
                },
                showInLegend: true
            }
        },
        series: [{
            name: 'Frequency',
            colorByPoint: true,
            data: plotData
        }]
    });
}

let pieSetValues = [];
let piePlotRepeats = [];

let pieRepeatSelect = new SlimSelect({
    select: "#pie-repeats-sel",
    placeholder: 'Select Repeats',
});

pieRepeatSelect.beforeOnChange = function(info) {
    if (info.length != 0) {
        let allValues = _.map(info, 'value')
        pieSetValues = allValues;
        let lastValue = info[info.length - 1].value
        if (lastValue == 'select-all') {
            for (let k in repeatSet) {
                pieSetValues = pieSetValues.concat([`select-all-${k}`])
                pieSetValues = pieSetValues.concat(repeatSet[k]);
                pieSetValues = _.uniq(pieSetValues);
                pieRepeatSelect.set([]);
                pieRepeatSelect.set(pieSetValues);
            }
        } else if (lastValue.slice(0, 10) == 'select-all') {
            let kmer = lastValue.slice(11, lastValue.length)
            pieSetValues = allValues.concat(repeatSet[kmer]);
            pieSetValues = _.uniq(pieSetValues);
            pieRepeatSelect.set([]);
            pieRepeatSelect.set(pieSetValues);
        }
        piePlotRepeats = _.filter(pieSetValues, function(d) {
            return d.slice(0, 10) != 'select-all';
        });
        piePlot(piePlotRepeats)
    }
}

pieRepeatSelect.onChange = function(info) {
    let currentValues = _.map(info, 'value');
    if (currentValues.length == 0 && pieSetValues.length - currentValues.length == 1 && pieSetValues[0].slice(0, 10) == 'select-all') {
        //pass
    } else if (pieSetValues.length - currentValues.length == 1) {
        let removedValue = (_.difference(pieSetValues, currentValues))[0];
        if (removedValue == 'select-all') {
            currentValues = [];
            pieRepeatSelect.set([]);
            pieRepeatSelect.set(currentValues);
            pieSetValues = currentValues;
            piePlotRepeats = _.filter(pieSetValues, function(d) {
                return d.slice(0, 10) != 'select-all';
            });
            piePlot('repeat', piePlotRepeats);
        } else if (removedValue.length > 10 && removedValue.slice(0, 10) == 'select-all') {
            let kmer = removedValue.slice(11, removedValue.length);
            let tempCurrentValues = _.difference(currentValues, repeatSet[kmer]);
            currentValues = tempCurrentValues;
            pieSetValues = currentValues;
            if (_.intersection(repeatSet[kmer], piePlotRepeats).length > 0) {
                pieRepeatSelect.set([]);
                pieRepeatSelect.set(currentValues);
                piePlotRepeats = _.filter(pieSetValues, function(d) {
                    return d.slice(0, 10) != 'select-all';
                });
                piePlot(piePlotRepeats);
            }
        } else {
            pieSetValues = currentValues;
            piePlotRepeats = _.filter(pieSetValues, function(d) {
                return d.slice(0, 10) != 'select-all';
            });
            piePlot(piePlotRepeats);
        }
    } else if (currentValues.length == 0) {
        //pass;
    }
}

const dataTypeRadios = document.getElementsByName('data-type');
for (let i in dataTypeRadios) {
    let radio = dataTypeRadios[i];
    if (radio.tagName == 'INPUT') {
        radio.onchange = function() { piePlot(piePlotRepeats); }
    }
}

const lenGroupRadios = document.getElementsByName('len-group');
for (let i in lenGroupRadios) {
    let radio = lenGroupRadios[i];
    if (radio.tagName == 'INPUT') {
        radio.onchange = function() { piePlot(piePlotRepeats); }
    }
}