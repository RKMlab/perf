/*-- LINE PLOT ------------------------------------------------------------------------------*/

const lineDatum = function(data, xdata, repeats, minL, maxL) {
    let datum = [];
    if (xdata == 0) {
        for (let r in repeats) {
            let repeat = repeats[r];
            datum.push({
                name: repeat,
                data: data['len'][repeat].slice(minL - minLength, maxL - minLength + 1)
            })
        }
    } else {
        for (let r in repeats) {
            let repeat = repeats[r];
            datum.push({
                name: repeat,
                data: data['unit'][repeat].slice(minL - minUnits, maxL - minUnits + 1)
            })
        }
    }
    return datum
}

let a = _.flatMap(plotInfo['len'], function(d) { return d.length + minLength - 1; });
document.getElementById('max-linex').value = _.max(a);

const linePlot = function(repeats) {
    const maxL = parseInt(document.getElementById('max-linex').value);
    const minL = parseInt(document.getElementById('min-linex').value);
    const xdata = document.getElementById('line-select').selectedIndex;
    let xLabel;
    if (xdata == 0) { xLabel = 'Sequence length(bp)'; } else { xLabel = 'Repeat Units'; }
    plotData = lineDatum(plotInfo, xdata, repeats, minL, maxL);

    Highcharts.chart('line-plot-svg', {
        chart: {
            marginTop: 40,
            marginLeft: 100
        },

        title: {
            text: null
        },

        yAxis: {
            title: {
                text: 'Frequency'
            }
        },

        xAxis: {
            title: {
                text: xLabel
            },
            labels: {
                format: '{value:.0f}'
            },
            tickInterval: 1
        },

        legend: {
            layout: 'vertical',
            align: 'right',
            verticalAlign: 'middle'
        },

        plotOptions: {
            series: {
                pointStart: parseInt(minL)
            }
        },

        series: plotData

    });
}

let lineSetValues = [];
let linePlotRepeats = [];

let lineRepeatSelect = new SlimSelect({
    select: "#line-repeats-sel",
    placeholder: 'Select Repeats',
});

lineRepeatSelect.beforeOnChange = function(info) {
    if (info.length != 0) {
        let allValues = _.map(info, 'value')
        lineSetValues = allValues;
        let lastValue = info[info.length - 1].value
        if (lastValue == 'select-all') {
            for (let k in repeatSet) {
                lineSetValues = lineSetValues.concat([`select-all-${k}`])
                lineSetValues = lineSetValues.concat(repeatSet[k]);
                lineSetValues = _.uniq(lineSetValues);
                lineRepeatSelect.set([]);
                lineRepeatSelect.set(lineSetValues);
            }
        } else if (lastValue.length > 10 && lastValue.slice(0, 10) == 'select-all') {
            let kmer = lastValue.slice(11, lastValue.length)
            lineSetValues = allValues.concat(repeatSet[kmer]);
            lineSetValues = _.uniq(lineSetValues);
            lineRepeatSelect.set([]);
            lineRepeatSelect.set(lineSetValues);
        }
        linePlotRepeats = _.filter(lineSetValues, function(d) {
            return d.slice(0, 10) != 'select-all';
        });
        linePlot(linePlotRepeats)
    }
}

lineRepeatSelect.onChange = function(info) {
    let currentValues = _.map(info, 'value');
    if (currentValues.length == 0 && lineSetValues.length - currentValues.length == 1 && lineSetValues[0].slice(0, 10) == 'select-all') {
        //pass
    } else if (lineSetValues.length - currentValues.length == 1) {
        let removedValue = (_.difference(lineSetValues, currentValues))[0];
        if (removedValue == 'select-all') {
            currentValues = [];
            lineRepeatSelect.set([]);
            lineRepeatSelect.set(currentValues);
            lineSetValues = currentValues;
            linePlotRepeats = [];
            linePlot(linePlotRepeats);
        } else if (removedValue.length > 10 && removedValue.slice(0, 10) == 'select-all') {
            let kmer = removedValue.slice(11, removedValue.length);
            let tempCurrentValues = _.difference(currentValues, repeatSet[kmer]);
            currentValues = tempCurrentValues;
            lineSetValues = currentValues;
            if (_.intersection(repeatSet[kmer], linePlotRepeats).length > 0) {
                lineRepeatSelect.set([]);
                lineRepeatSelect.set(currentValues);
                linePlotRepeats = _.filter(lineSetValues, function(d) {
                    return d.slice(0, 10) != 'select-all';
                });
                linePlot(linePlotRepeats);
            }
        } else {
            lineSetValues = currentValues;
            linePlotRepeats = _.filter(lineSetValues, function(d) {
                return d.slice(0, 10) != 'select-all';
            });
            linePlot(linePlotRepeats);
        }
    } else if (currentValues.length == 0) {
        //pass;
    }
}

document.getElementById('min-linex').oninput = function() { linePlot(linePlotRepeats) };
document.getElementById('max-linex').oninput = function() { linePlot(linePlotRepeats) };
document.getElementById('line-select').onchange = function() {
    if (this.selectedIndex == 0) {
        document.getElementById('min-linex').value = minLength
    } else {
        document.getElementById('min-linex').value = minUnits
    }
    linePlot(linePlotRepeats);
};