/*-- LINE PLOT ------------------------------------------------------------------------------*/

const lineDatum = function(data, repeats, minL, maxL) {
    let datum = [];
    for (let r in repeats) {
        let repeat = repeats[r];
        datum.push({
            name: repeat,
            data: data[repeat].slice(minL - 12, maxL - 12 + 1)
        })
    }
    return datum
}

let a = _.flatMap(plotInfo, function(d) { return d.length + 12 - 1; });
document.getElementById('max-length').value = _.max(a);

const linePlot = function(repeats) {
    const maxL = document.getElementById('max-length').value;
    const minL = document.getElementById('min-length').value;
    console.log()
    let plotData = lineDatum(plotInfo, repeats, minL, maxL);

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
                text: 'Sequence length(bp)'
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
        let allValues = _.map(info, function(d) {
            return d.value;
        })
        lineSetValues = allValues;
        let lastValue = info[info.length - 1].value
        if (lastValue == 'select-all') {
            for (let k in repeatSet) {
                lineSetValues = lineSetValues.concat(repeatSet[k]);
                lineSetValues = _.uniq(lineSetValues);
                lineRepeatSelect.set([]);
                lineRepeatSelect.set(lineSetValues);
            }
        } else if (lastValue.slice(0, 10) == 'select-all') {
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
    let currentValues = _.map(info, function(d) { return d.value; })
    if (lineSetValues.length - currentValues.length == 1) {
        let removedValue = (_.difference(lineSetValues, currentValues))[0];
        if (removedValue == 'select-all') {
            currentValues = [];
        } else if (removedValue.slice(0, 10) == 'select-all') {
            let kmer = removedValue.slice(11, removedValue.length)
            let tempCurrentValues = _.difference(currentValues, repeatSet[kmer]);
            currentValues = tempCurrentValues;
        }
        lineRepeatSelect.set([])
        lineRepeatSelect.set(currentValues)
        lineSetValues = currentValues;
        linePlotRepeats = _.filter(lineSetValues, function(d) {
            return d.slice(0, 10) != 'select-all';
        });
        linePlot(linePlotRepeats)
    }
}
document.getElementById('min-length').oninput = function() { linePlot(linePlotRepeats) };
document.getElementById('max-length').oninput = function() { linePlot(linePlotRepeats) };