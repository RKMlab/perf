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
    const header = ['Repeat Class', 'Frequency', 'Frequency per million', '% Frequency', 'Reads', 'Reads per million', '% Reads', 'Bases']
    header.forEach(function(e){ const headCell = document.createElement('th'); headCell.innerHTML = e; tableHeadRow.appendChild(headCell); })
    tableHead.appendChild(tableHeadRow);

    const tableBody = document.createElement('tbody');

    // const totalRepBases = _.sum(_.map(Object.keys(tableData), o => { return tableData[o][0]; }));
    const totalRepFreq = tableData.totalRepFreq;
    const totalRepReads = tableData.totalRepReads;
    const totals = [totalRepReads, totalRepFreq]
    const cell_keys = ['instances']
    for (let rep in allRepClasses) {
        rep = allRepClasses[rep];
        const row = document.createElement('tr');
        const rep_cell = document.createElement('td');
        rep_cell.innerHTML = rep; row.appendChild(rep_cell);

        let rowData = [];

        if (rep in tableData.repClassInfo) {
            const repInfo = tableData.repClassInfo[rep];
            rowData.push(repInfo['instances']);
            rowData.push(repInfo['instances_norm']);
            rowData.push(((repInfo['instances']/totalRepFreq)*100).toFixed(3));

            rowData.push(repInfo['reads']);
            rowData.push(repInfo['reads_norm']);
            rowData.push(((repInfo['reads']/totalRepReads)*100).toFixed(3));

            rowData.push(repInfo['bases']);
        }
        else { rowData = Array(7).fill(0); }
        rowData.forEach(function(e){ const cell = document.createElement('td'); cell.innerHTML = e; row.appendChild(cell); })

        tableBody.appendChild(row);
    }

    table.appendChild(tableHead);
    table.appendChild(tableBody);
    tableDOM.appendChild(table);
}

updateSummaryTableData('rep-summary-table', data.info.repInfo);