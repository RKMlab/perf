/*
    The "tables.js" contains JS code which updates the data for
    tables present in the PERF-analysis module.

    This web application is developed with Semantic-ui frame work.
    The charts are build using Apex-Chart js charting library.

    All the data for the report is derived from analyse_data.js
    data = {info: {genomeInfo: {}, repInfo: {}, plotInfo: {}}}

*/



const updateSummaryTableData = function(tableId, tableData) {
    const tableDOM = document.getElementById(tableId);

    const table = document.createElement('table');
    table.className = "ui sortable celled table";
    const tableHead = document.createElement('thead')
    const tableHeadRow = document.createElement('tr');
    const header = ['Repeat Class', 'Frequency', '% Frequency', 'Bases', '% Bases']
    header.forEach(function(e){ const headCell = document.createElement('th'); headCell.innerHTML = e; tableHeadRow.appendChild(headCell); })
    tableHead.appendChild(tableHeadRow);

    const tableBody = document.createElement('tbody');
    const totalRepBases = _.sum(_.map(Object.keys(tableData), o => { return tableData[o][0]; }));
    const totalRepFreq = _.sum(_.map(Object.keys(tableData), o => { return tableData[o][1]; }));
    const totals = [totalRepBases, totalRepFreq]
    for (let rep in allRepClasses) {
        rep = allRepClasses[rep];
        const row = document.createElement('tr');
        const rep_cell = document.createElement('td');
        rep_cell.innerHTML = rep; row.appendChild(rep_cell);

        const rowData = []; 
        tableData[rep].forEach(function(d, i){ rowData.push(d); rowData.push(((d/totals[i])*100).toFixed(3)); });
        rowData.forEach(function(e){ const cell = document.createElement('td'); cell.innerHTML = e; row.appendChild(cell); })

        tableBody.appendChild(row);
    }

    table.appendChild(tableHead);
    table.appendChild(tableBody);
    tableDOM.appendChild(table);
}

updateSummaryTableData('rep-summary-table', barData);

const updateLongestRepeatsTableData = function(tableId, tableData) {
    const tableDOM = document.getElementById(tableId);
    const table = document.createElement('table');
    table.className = "ui sortable celled table";
    const tableHead = document.createElement('thead')
    const tableHeadRow = document.createElement('tr');
    const header = ['Sequence id', 'Start', 'Stop', 'Repeat Class', 'Repeat length', 'Strand', 'Units', 'Actual Repeat'];
    header.forEach(function(e){ const headCell = document.createElement('th'); headCell.innerHTML = e; tableHeadRow.appendChild(headCell); })
    tableHead.appendChild(tableHeadRow);

    const tableBody = document.createElement('tbody');
    for (let d in tableData) {
        d = tableData[d];
        const row = document.createElement('tr');
        const rowData = ["seq", "start", "end", "repClass", "repLength", "repOri", "repUnit", "actualRep"];
        rowData.forEach(function(e){ const cell = document.createElement('td'); cell.innerHTML = d[e]; row.appendChild(cell); })
        tableBody.appendChild(row);
    }

    table.appendChild(tableHead);
    table.appendChild(tableBody);
    tableDOM.appendChild(table);
}

updateLongestRepeatsTableData('longest-repeats-table', data.info.repInfo.longestRepeats);
updateLongestRepeatsTableData('mostunits-repeats-table', data.info.repInfo.mostRepeatUnits);