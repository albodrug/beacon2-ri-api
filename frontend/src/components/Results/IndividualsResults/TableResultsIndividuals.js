import './TableResultsIndividuals.css'
import * as React from 'react'
import { useState, useEffect } from 'react'
import CrossQueries from '../../CrossQueries/CrossQueries'
import { FaBars, FaEye, FaEyeSlash } from 'react-icons/fa' // Import icons from react-icons library
import { FiLayers, FiDownload } from 'react-icons/fi'

function TableResultsIndividuals (props) {
  const [showDatsets, setShowDatasets] = useState(false)
  const [showResults, setShowResults] = useState(false)
  const [resultsSelected, setResultsSelected] = useState(props.results)
  const [arrayBeaconsIds, setArrayBeaconsIds] = useState([])
  const [errorMessage, setErrorMessage] = useState('')
  const [resultsSelectedFinal, setResultsSelectedFinal] = useState([])
  const [editable, setEditable] = useState([])
  const [trigger, setTrigger] = useState(false)
  const [trigger2, setTrigger2] = useState(false)
  const [exportMenuVisible, setExportMenuVisible] = useState(false)
  const [showCrossQuery, setShowCrossQuery] = useState(false)
  const [parameterCrossQuery, setParamCrossQuery] = useState('')

  const [note, setNote] = useState('')
  const [isOpenModal2, setIsOpenModal2] = useState(false)

  const [filterValues, setFilterValues] = useState({
    IndividualId: '',
    ethnicity: '',
    Beacon: '',
    interventionsOrProcedures: '',
    sex: '',
    diseases: '',
    treatments: '',
    phenotypicFeatures: ''
    // Add other column names here
  })

  const [menuVisible, setMenuVisible] = useState(false)

  const toggleMenu = () => {
    setMenuVisible(prevState => !prevState)
  }

  const [columnVisibility, setColumnVisibility] = useState({
    IndividualId: true,
    ethnicity: true,
    Beacon: true,
    interventionsOrProcedures: true,
    sex: true,
    diseases: true,
    treatments: true,
    phenotypicFeatures: true
    // Add more columns as needed
  })
  const showAllColumns = () => {
    const columns = document.querySelectorAll('th')
    const rows = document.querySelectorAll('td')

    // Update column visibility state and remove hidden class for all columns
    columns.forEach(column => {
      column.classList.remove('hidden')
      const columnName = column.dataset.columnName
      setColumnVisibility(prevState => ({
        ...prevState,
        [columnName]: true
      }))
    })

    // Change the icon of all rows to the normal eye
    rows.forEach(row => {
      row.classList.remove('hidden')
    })

    setColumnVisibility(prevState => {
      const updatedVisibility = {}
      Object.keys(prevState).forEach(column => {
        updatedVisibility[column] = true
      })
      return updatedVisibility
    })
  }

  const [filteredData, setFilteredData] = useState(editable)
  const toggleColumnVisibility = columnName => {
    const columns = document.querySelectorAll('th[data-column-name]')
    const rows = document.querySelectorAll(
      `td[data-column-name="${columnName}"]`
    )

    columns.forEach(column => {
      if (column.dataset.columnName === columnName) {
        column.classList.toggle('hidden')
      }
    })

    rows.forEach(row => {
      row.classList.toggle('hidden')
    })

    setColumnVisibility(prevVisibility => ({
      ...prevVisibility,
      [columnName]: !prevVisibility[columnName]
    }))
  }

  const handleFilterChange = (e, columnName) => {
    const { value } = e.target
    setFilterValues({ ...filterValues, [columnName]: value })

    const updatedFilteredData = editable.filter(row =>
      row[columnName].toLowerCase().includes(value.toLowerCase())
    )

    setFilteredData(updatedFilteredData)
  }

  const toggleExportMenu = () => {
    setExportMenuVisible(prevState => !prevState)
  }

  const exportToCSV = () => {
    // Ensure props.results is not null or undefined
    if (!props.results) return;
  
    // Get all keys from the first row of props.results
    const header = Object.keys(props.results[0]);
  
    // Convert each row to CSV format
    const csv = [
      header.join(','), // Header row
      ...props.results.map(row =>
        header.map(fieldName => {
          const value = row[fieldName];
          // Check if the value is an object
          if (typeof value === 'object') {
            // Stringify the object
            return JSON.stringify(value);
          } else {
            // Otherwise, return the value as is
            return value;
          }
        }).join(',')
      )
    ].join('\n');
  
    // Create a blob object from the CSV content
    const blob = new Blob([csv], { type: 'text/csv' });
  
    // Create a URL for the blob object
    const url = window.URL.createObjectURL(blob);
  
    // Create a temporary <a> element to trigger the download
    const link = document.createElement('a');
    link.href = url;
    link.setAttribute('download', 'exported_data.csv');
  
    // Programmatically click the link to start the download
    document.body.appendChild(link);
    link.click();
  
    // Clean up by revoking the URL and removing the temporary <a> element
    window.URL.revokeObjectURL(url);
    document.body.removeChild(link);
  };
  

  const exportToJSON = () => {
    // Ensure props.results is not null or undefined
    if (!props.results) return;
  
    // Convert the results to JSON
    const jsonString = JSON.stringify(props.results, null, 2);
  
    // Create a blob object from the JSON content
    const blob = new Blob([jsonString], { type: 'application/json' });
  
    // Create a URL for the blob object
    const url = URL.createObjectURL(blob);
  
    // Create a temporary <a> element to trigger the download
    const link = document.createElement('a');
    link.href = url;
    link.setAttribute('download', 'exported_data.json');
  
    // Programmatically click the link to start the download
    document.body.appendChild(link);
    link.click();
  
    // Clean up by revoking the URL and removing the temporary <a> element
    URL.revokeObjectURL(url);
    document.body.removeChild(link);
  };
  


  const showNote = e => {
    setNote(e)
    setIsOpenModal2(true)
  }

  const handleShowCrossQuery = e => {
    setShowCrossQuery(true)
    console.log(e.target.innerText)
    setParamCrossQuery(e.target.innerText)
  }

  useEffect(() => {
    if (props.show === 'full') {
      setResultsSelectedFinal(resultsSelected)
      setShowResults(true)
      setShowDatasets(false)
      setTrigger(true)
    }

    if (resultsSelected.length === 0) {
      setErrorMessage('NO RESULTS')
    }
    resultsSelected.forEach((element, index) => {
      arrayBeaconsIds.push(element[0])
    })

    resultsSelectedFinal.forEach((element, index) => {
      if (element[1] !== undefined) {
        let eth_id = ''
        let eth_label = ''
        let stringEth = ''

        if (element[1].ethnicity !== '' && element[1].ethnicity !== undefined) {
          if (element[1].ethnicity.id !== undefined) {
            eth_id = element[1].ethnicity.id
          }

          eth_label = element[1].ethnicity.label
          stringEth = `${eth_id} / ${eth_label} `
        } else {
          stringEth = ''
        }

        let sex_id = ''
        let sex_label = ''
        let stringSex = ''

        if (element[1].sex !== '') {
          sex_id = element[1].sex.id
          sex_label = element[1].sex.label
          stringSex = `${element[1].sex.label} / ${element[1].sex.id}`
        } else {
          stringSex = ''
        }

        let geographic_id = ''
        let geographic_label = ''
        let stringGeographic = ''

        if (
          element[1].geographicOrigin !== '' &&
          element[1].geographicOrigin !== undefined
        ) {
          geographic_id = element[1].geographicOrigin.id
          geographic_label = element[1].geographicOrigin.label
          stringGeographic = `${geographic_id} / ${geographic_label}`
        } else {
          stringGeographic = ''
        }

        let measuresJson = ''
        if (element[1].measures !== '' && element[1].measures !== undefined) {
          if (typeof element[1].measures === 'object') {
            measuresJson = JSON.stringify(element[1].measures, null, 2)
              .replaceAll('[', '')
              .replaceAll(']', '')
              .replaceAll('{', '')
              .replaceAll('}', '')
              .replaceAll(',', '')
              .replaceAll(' ,', '')
              .replaceAll(', ', '')
              .replaceAll('"', '')
          } else {
            measuresJson = JSON.stringify(element[1].measures, null, 2)
              .replaceAll('[', '')
              .replaceAll(']', '')
              .replaceAll('{', '')
              .replaceAll('}', '')
              .replaceAll(',', '')
              .replaceAll(' ,', '')
              .replaceAll(', ', '')
              .replaceAll('"', '')
          }
        }

        let phenoJson = ''
        if (
          element[1].phenotypicFeatures !== '' &&
          element[1].phenotypicFeatures !== undefined
        ) {
          if (typeof element[1].phenotypicFeatures === 'object') {
            phenoJson = JSON.stringify(element[1].phenotypicFeatures, null, 2)
              .replaceAll('[', '')
              .replaceAll(']', '')
              .replaceAll('{', '')
              .replaceAll('}', '')
              .replaceAll(',', '')
              .replaceAll(' ,', '')
              .replaceAll(', ', '')
              .replaceAll('"', '')
          } else {
            phenoJson = JSON.stringify(element[1].phenotypicFeatures, null, 2)
              .replaceAll('[', '')
              .replaceAll(']', '')
              .replaceAll('{', '')
              .replaceAll('}', '')
              .replaceAll(',', '')
              .replaceAll(' ,', '')
              .replaceAll(', ', '')
              .replaceAll('"', '')
          }
        }

        let interventionsProcedures = ''

        if (
          element[1].interventionsOrProcedures !== '' &&
          element[1].interventionsOrProcedures !== undefined
        ) {
          if (typeof element[1].interventionsOrProcedures === 'object') {
            interventionsProcedures = JSON.stringify(
              element[1].interventionsOrProcedures,
              null,
              2
            )
              .replaceAll('[', '')
              .replaceAll(']', '')
              .replaceAll('{', '')
              .replaceAll('}', '')
              .replaceAll(',', '')
              .replaceAll(' ,', '')
              .replaceAll(', ', '')
              .replaceAll('"', '')
          } else {
            interventionsProcedures = JSON.stringify(
              element[1].interventionsOrProcedures,
              null,
              2
            )
              .replaceAll('[', '')
              .replaceAll(']', '')
              .replaceAll('{', '')
              .replaceAll('}', '')
              .replaceAll(',', '')
              .replaceAll(' ,', '')
              .replaceAll(', ', '')
              .replaceAll('"', '')
          }
        }

        let diseases = ''

        if (element[1].diseases !== '' && element[1].diseases !== undefined) {
          if (typeof element[1].diseases === 'object') {
            diseases = JSON.stringify(element[1].diseases, null, 2)
              .replaceAll('[', '')
              .replaceAll(']', '')
              .replaceAll('{', '')
              .replaceAll('}', '')
              .replaceAll(',', '')
              .replaceAll(' ,', '')
              .replaceAll(', ', '')
              .replaceAll('"', '')
          } else {
            diseases = JSON.stringify(element[1].diseases, null, 2)
              .replaceAll('[', '')
              .replaceAll(']', '')
              .replaceAll('{', '')
              .replaceAll('}', '')
              .replaceAll(',', '')
              .replaceAll(' ,', '')
              .replaceAll(', ', '')
              .replaceAll('"', '')
          }
        }

        let treatments = ''

        if (
          element[1].treatments !== '' &&
          element[1].treatments !== undefined
        ) {
          if (typeof element[1].treatments === 'object') {
            treatments = JSON.stringify(element[1].treatments, null, 2)
              .replaceAll('[', '')
              .replaceAll(']', '')
              .replaceAll('{', '')
              .replaceAll('}', '')
              .replaceAll(',', '')
              .replaceAll(' ,', '')
              .replaceAll(', ', '')
              .replaceAll('"', '')
          } else {
            treatments = JSON.stringify(element[1].treatments, null, 2)
              .replaceAll('[', '')
              .replaceAll(']', '')
              .replaceAll('{', '')
              .replaceAll('}', '')
              .replaceAll(',', '')
              .replaceAll(' ,', '')
              .replaceAll(', ', '')
              .replaceAll('"', '')
          }
        }

        editable.push({
          id: index,
          IndividualId: element[1].id,
          Beacon: element[0],
          ethnicity: stringEth,
          interventionsOrProcedures: interventionsProcedures,
          sex: stringSex,
          diseases: diseases,
          treatments: treatments,
          phenotypicFeatures: phenoJson
        })

        if (index === resultsSelectedFinal.length - 1) {
          setTrigger2(true)
        }
      }
    })
  }, [trigger, resultsSelectedFinal])

  useEffect(() => {
    setShowDatasets(true)
  }, [])

  return (
    <div className='containerBeaconResults'>
      {showDatsets === true &&
        props.beaconsList.map(result => {
          return (
            <>
              {props.show !== 'full' &&
                props.resultsPerDataset.map((element, index) => {
                  return (
                    <>
                      {element[1][index] === true &&
                        props.show === 'boolean' && (
                          <h6 className='foundResult'>YES</h6>
                        )}
                      {element[1][index] === false &&
                        props.show === 'boolean' && (
                          <h5 className='NotFoundResult'>No, sorry</h5>
                        )}
                      {props.show === 'count' &&
                        element[2][index] !== 0 &&
                        element[2][index] !== 1 && (
                          <h6 className='foundResult'>
                            {element[2][index]} RESULTS
                          </h6>
                        )}
                      {props.show === 'count' && element[2][index] === 0 && (
                        <h5 className='NotFoundResult'>
                          {element[2][index]} RESULTS
                        </h5>
                      )}
                      {props.show === 'count' && element[2][index] === 1 && (
                        <h6 className='foundResult'>
                          {element[2][index]} RESULT
                        </h6>
                      )}
                    </>
                  )
                })}
            </>
          )
        })}

      {!showCrossQuery &&
        showDatsets === false &&
        showResults === true &&
        trigger2 === true && (
          <div className='table-container'>
            <div className='menu-icon-container'>
              <div className='export-menu'>
                <button className='exportButton' onClick={toggleExportMenu}>
                  <FiDownload />
                </button>
                {exportMenuVisible && (
                  <>
                    <ul className='column-list'>
                      <li onClick={exportToJSON}>Export to JSON</li>
                      <li onClick={exportToCSV}>Export to CSV</li>
                    </ul>
                  </>
                )}
              </div>
              <div className='menu-container'>
                <FaBars onClick={toggleMenu} />
                {menuVisible && (
                  <>
                    <ul className='column-list'>
                      <li onClick={showAllColumns}>
                        Show All Columns
                        <FiLayers />
                      </li>
                      {Object.keys(columnVisibility).map(column => (
                        <li
                          key={column}
                          onClick={() => toggleColumnVisibility(column)}
                        >
                          {column}
                          {columnVisibility[column] ? (
                            <FaEye />
                          ) : (
                            <FaEyeSlash />
                          )}
                        </li>
                      ))}
                    </ul>
                  </>
                )}
              </div>
            </div>
            <div className='header-container'>
              <table className='tableResults'>
                <thead className='theadResults'>
                  <tr>
                    <th
                      className={`sticky-header ${
                        columnVisibility.IndividualId ? 'visible' : 'hidden'
                      }`}
                    >
                      <span>Individual ID</span>
                      <button
                        onClick={() => toggleColumnVisibility('IndividualId')}
                      >
                        {columnVisibility.IndividualId ? (
                          <FaEye />
                        ) : (
                          <FaEyeSlash />
                        )}
                      </button>
                      <input
                        type='text'
                        placeholder='Filter Individual ID'
                        onChange={e => handleFilterChange(e, 'IndividualId')}
                      />
                    </th>
                    <th
                      className={`sticky-header ${
                        columnVisibility.ethnicity ? 'visible' : 'hidden'
                      }`}
                    >
                      <span>Ethnicity</span>
                      <button
                        onClick={() => toggleColumnVisibility('ethnicity')}
                      >
                        {columnVisibility.ethnicity ? (
                          <FaEye />
                        ) : (
                          <FaEyeSlash />
                        )}
                      </button>
                      <input
                        type='text'
                        placeholder='Filter Ethnicity'
                        onChange={e => handleFilterChange(e, 'ethnicity')}
                      />
                    </th>
                    <th
                      className={`sticky-header ${
                        columnVisibility.Beacon ? 'visible' : 'hidden'
                      }`}
                    >
                      <span>Beacon</span>
                      <button onClick={() => toggleColumnVisibility('Beacon')}>
                        {columnVisibility.Beacon ? <FaEye /> : <FaEyeSlash />}
                      </button>
                      <input
                        type='text'
                        placeholder='Filter Beacon'
                        onChange={e => handleFilterChange(e, 'Beacon')}
                      />
                    </th>
                    <th
                      className={`sticky-header ${
                        columnVisibility.interventionsOrProcedures
                          ? 'visible'
                          : 'hidden'
                      }`}
                    >
                      <span>Procedures</span>
                      <button
                        onClick={() =>
                          toggleColumnVisibility('interventionsOrProcedures')
                        }
                      >
                        {columnVisibility.interventionsOrProcedures ? (
                          <FaEye />
                        ) : (
                          <FaEyeSlash />
                        )}
                      </button>
                      <input
                        type='text'
                        placeholder='Filter Procedures'
                        onChange={e =>
                          handleFilterChange(e, 'interventionsOrProcedures')
                        }
                      />
                    </th>
                    <th
                      className={`sticky-header ${
                        columnVisibility.sex ? 'visible' : 'hidden'
                      }`}
                    >
                      <span>Sex</span>
                      <button onClick={() => toggleColumnVisibility('sex')}>
                        {columnVisibility.sex ? <FaEye /> : <FaEyeSlash />}
                      </button>
                      <input
                        type='text'
                        placeholder='Filter Sex'
                        onChange={e => handleFilterChange(e, 'sex')}
                      />
                    </th>
                    <th
                      className={`sticky-header ${
                        columnVisibility.diseases ? 'visible' : 'hidden'
                      }`}
                    >
                      <span>Diseases</span>
                      <button
                        onClick={() => toggleColumnVisibility('diseases')}
                      >
                        {columnVisibility.diseases ? <FaEye /> : <FaEyeSlash />}
                      </button>
                      <input
                        type='text'
                        placeholder='Filter Diseases'
                        onChange={e => handleFilterChange(e, 'diseases')}
                      />
                    </th>
                    <th
                      className={`sticky-header ${
                        columnVisibility.treatments ? 'visible' : 'hidden'
                      }`}
                    >
                      <span>Treatments</span>
                      <button
                        onClick={() => toggleColumnVisibility('treatments')}
                      >
                        {columnVisibility.treatments ? (
                          <FaEye />
                        ) : (
                          <FaEyeSlash />
                        )}
                      </button>
                      <input
                        type='text'
                        placeholder='Filter Treatments'
                        onChange={e => handleFilterChange(e, 'treatments')}
                      />
                    </th>
                    <th
                      className={`sticky-header ${
                        columnVisibility.phenotypicFeatures
                          ? 'visible'
                          : 'hidden'
                      }`}
                    >
                      <span>Phenotypic Features</span>
                      <button
                        onClick={() =>
                          toggleColumnVisibility('phenotypicFeatures')
                        }
                      >
                        {columnVisibility.phenotypicFeatures ? (
                          <FaEye />
                        ) : (
                          <FaEyeSlash />
                        )}
                      </button>
                      <input
                        type='text'
                        placeholder='Filter Phenotypic Features'
                        onChange={e =>
                          handleFilterChange(e, 'phenotypicFeatures')
                        }
                      />
                    </th>

                    {/* Add more column headers here */}
                  </tr>
                </thead>
              </table>
            </div>
            <div className='body-container'>
              <table className='tableResults'>
                <tbody>
                  {filteredData.map((row, index) => (
                    <tr key={index}>
                      <td
                        className={
                          columnVisibility.IndividualId ? 'visible' : 'hidden'
                        }
                      >
                        {row.IndividualId}
                      </td>
                      <td
                        className={
                          columnVisibility.ethnicity ? 'visible' : 'hidden'
                        }
                      >
                        {row.ethnicity}
                      </td>
                      <td
                        className={
                          columnVisibility.Beacon ? 'visible' : 'hidden'
                        }
                      >
                        {row.Beacon}
                      </td>
                      <td
                        className={
                          columnVisibility.interventionsOrProcedures
                            ? 'visible'
                            : 'hidden'
                        }
                      >
                        {row.interventionsOrProcedures}
                      </td>
                      <td
                        className={columnVisibility.sex ? 'visible' : 'hidden'}
                      >
                        {row.sex}
                      </td>
                      <td
                        className={
                          columnVisibility.diseases ? 'visible' : 'hidden'
                        }
                      >
                        {row.diseases}
                      </td>
                      <td
                        className={
                          columnVisibility.treatments ? 'visible' : 'hidden'
                        }
                      >
                        {row.treatments}
                      </td>
                      <td
                        className={
                          columnVisibility.phenotypicFeatures
                            ? 'visible'
                            : 'hidden'
                        }
                      >
                        {row.phenotypicFeatures}
                      </td>

                      {/* Render other row cells here */}
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </div>
        )}
      {showCrossQuery && (
        <CrossQueries
          parameter={parameterCrossQuery}
          collection={'individuals'}
          setShowCrossQuery={setShowCrossQuery}
        />
      )}
    </div>
  )
}

export default TableResultsIndividuals
