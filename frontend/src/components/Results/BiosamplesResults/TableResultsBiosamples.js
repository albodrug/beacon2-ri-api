import './TableResultsBiosamples.css'
import '../IndividualsResults/TableResultsIndividuals.css'
import '../../Dataset/BeaconInfo'
import * as React from 'react'
import {
  DataGrid,
  GridToolbar,
  selectedGridRowsSelector,
  gridFilteredSortedRowIdsSelector,
  GridToolbarContainer
} from '@mui/x-data-grid'
import { useState, useEffect } from 'react'
import CrossQueries from '../../CrossQueries/CrossQueries'

function CustomToolbar () {
  return <GridToolbarContainer></GridToolbarContainer>
}
function TableResultsBiosamples (props) {
  const [showDatsets, setShowDatasets] = useState(false)

  const [showResults, setShowResults] = useState(false)

  const [arrayBeaconsIds, setArrayBeaconsIds] = useState([])
  const [rows, setRows] = useState([])
  const [ids, setIds] = useState([])

  const [showCrossQuery, setShowCrossQuery] = useState(false)
  const [parameterCrossQuery, setParamCrossQuery] = useState('')

  const [errorMessage, setErrorMessage] = useState('')

  const [beaconsArrayResults, setBeaconsArrayResults] = useState([])

  const [beaconsArrayResultsOrdered, setBeaconsArrayResultsOrdered] = useState(
    []
  )

  const [resultsSelected, setResultsSelected] = useState(props.results)
  const [resultsSelectedFinal, setResultsSelectedFinal] = useState([])

  const [openDatasetArray, setOpenDataset] = useState([])
  const [openDatasetArray2, setOpenDataset2] = useState([])

  const [editable, setEditable] = useState([])

  const [trigger, setTrigger] = useState(false)
  const [trigger2, setTrigger2] = useState(false)

  const [triggerArray, setTriggerArray] = useState([])
  const [triggerArray2, setTriggerArray2] = useState([])

  const getSelectedRowsToExport = ({ apiRef }) => {
    const selectedRowIds = selectedGridRowsSelector(apiRef)
    if (selectedRowIds.size > 0) {
      return Array.from(selectedRowIds.keys())
    }

    return gridFilteredSortedRowIdsSelector(apiRef)
  }

  const handleShowCrossQuery = e => {
    setShowCrossQuery(true)
    console.log(e.target.innerText)
    setParamCrossQuery(e.target.innerText)
  }

  let columns = [
    {
      field: 'BiosampleId',
      headerName: 'Biosample ID',
      width: 150,
      headerClassName: 'super-app-theme--header',
      renderCell: params => (
        <button className='buttonId' onClick={handleShowCrossQuery}>
          {params.row.BiosampleId}
        </button>
      )
    },
    {
      field: 'Beacon',
      headerName: 'Beacon ID',
      width: 340,
      headerClassName: 'super-app-theme--header'
    },
    {
      field: 'individualId',
      headerName: 'Individual ID',
      width: 150,
      headerClassName: 'super-app-theme--header'
    },
    {
      field: 'biosampleStatus',
      headerName: 'Biosample status',
      width: 240,
      headerClassName: 'super-app-theme--header'
    },
    // {
    //   field: 'collectionDate',
    //   headerName: 'Collection date',
    //   width: 250,
    //   headerClassName: 'super-app-theme--header'
    // },
    // {
    //   field: 'collectionMoment',
    //   headerName: 'Collection moment',
    //   width: 350,
    //   headerClassName: 'super-app-theme--header'
    // },
    {
      field: 'sampleOriginType',
      headerName: 'Sample origin type',
      width: 350,
      headerClassName: 'super-app-theme--header',
      cellClass: 'pre'
    },
    // {
    //   field: 'sampleOriginDetail',
    //   headerName: 'Sample origin detail',
    //   width: 200,
    //   headerClassName: 'super-app-theme--header'
    // },
    // {
    //   field: 'obtentionProcedure',
    //   headerName: 'Obtention procedure',
    //   width: 300,
    //   headerClassName: 'super-app-theme--header'
    // },
    // {
    //   field: 'tumorProgression',
    //   headerName: 'Tumor progression',
    //   width: 350,
    //   headerClassName: 'super-app-theme--header'
    // },
    // {
    //   field: 'tumorGrade',
    //   headerName: 'Tumor Grade',
    //   width: 200,
    //   headerClassName: 'super-app-theme--header'
    // },
    // {
    //   field: 'pathologicalStage',
    //   headerName: 'Pathological stage',
    //   width: 350,
    //   headerClassName: 'super-app-theme--header'
    // },
    {
      field: 'pathologicalTnmFinding',
      headerName: 'Pathological TNM findings',
      width: 300,
      headerClassName: 'super-app-theme--header'
    },
    {
      field: 'histologicalDiagnosis',
      headerName: 'Histological diagnosis',
      width: 350,
      headerClassName: 'super-app-theme--header'
    },
    // {
    //   field: 'diagnosticMarkers',
    //   headerName: 'Diagnostic markers',
    //   width: 300,
    //   headerClassName: 'super-app-theme--header'
    // },
    // {
    //   field: 'phenotypicFeatures',
    //   headerName: 'Phenotypic features',
    //   width: 300,
    //   headerClassName: 'super-app-theme--header'
    // },
    {
      field: 'measurements',
      headerName: 'Measurements',
      width: 300,
      headerClassName: 'super-app-theme--header'
    },
    // {
    //   field: 'sampleProcessing',
    //   headerName: 'Sample processing',
    //   width: 300,
    //   headerClassName: 'super-app-theme--header'
    // },
    {
      field: 'sampleStorage',
      headerName: 'Sample storage',
      width: 300,
      headerClassName: 'super-app-theme--header'
    }
  ]

  const handleSeeResults = e => {
    resultsSelected.forEach(element => {
      if (element[0] === e) {
        resultsSelectedFinal.push(element)
      }
    })
    setShowResults(true)
    setShowDatasets(false)
    setTrigger(true)
  }

  function getOccurrence (array, value) {
    var count = 0
    array.forEach(v => v === value && count++)
    return count
  }

  useEffect(() => {
    if (props.show === 'full') {
      setResultsSelectedFinal(resultsSelected)
      setShowResults(true)
      setShowDatasets(false)
      setTrigger(true)
    }
    setRows([])
    setIds([])
    if (resultsSelected.length === 0) {
      setErrorMessage('NO RESULTS')
    }
    resultsSelected.forEach((element, index) => {
      arrayBeaconsIds.push(element[0])
    })
    resultsSelectedFinal.forEach((element, index) => {
      if (element[1] !== undefined) {
        let biosampleStatus_id = ''
        let biosampleStatus_label = ''
        let stringBiosampleStatus = ''

        if (
          element[1].biosampleStatus !== '' &&
          element[1].biosampleStatus !== undefined
        ) {
          if (element[1].biosampleStatus.id !== undefined) {
            biosampleStatus_id = element[1].biosampleStatus.id
          }
          if (element[1].biosampleStatus.label !== undefined) {
            biosampleStatus_label = element[1].biosampleStatus.label
          }

          stringBiosampleStatus = `${biosampleStatus_id} / ${biosampleStatus_label} `
        } else {
          stringBiosampleStatus = ''
        }

        let sampleOriginType_id = ''
        let sampleOriginType_label = ''
        let stringSampleOriginType = ''

        if (
          element[1].sampleOriginType !== '' &&
          element[1].sampleOriginType !== undefined
        ) {
          sampleOriginType_id = element[1].sampleOriginType.id
          sampleOriginType_label = element[1].sampleOriginType.label
          stringSampleOriginType = `${element[1].sampleOriginType.label} / ${element[1].sampleOriginType.id}`
        } else {
          stringSampleOriginType = ''
        }

        let sampleOriginDetail_id = ''
        let sampleOriginDetail_label = ''
        let stringSampleOriginDetail = ''

        if (
          element[1].sampleOriginDetail !== '' &&
          element[1].sampleOriginDetail !== undefined
        ) {
          sampleOriginDetail_id = element[1].sampleOriginDetail.id
          sampleOriginDetail_label = element[1].sampleOriginDetail.label
          stringSampleOriginDetail = `${sampleOriginDetail_id} / ${sampleOriginDetail_label}`
        } else {
          stringSampleOriginDetail = ''
        }

        let collectionDateJson = []
        if (
          element[1].collectionDate !== '' &&
          element[1].collectionDate !== undefined
        ) {
          if (typeof element[1].collectionDate === 'string') {
            collectionDateJson = element[1].collectionDate
          } else {
            collectionDateJson = element[1].collectionDate.toString()
          }
        }

        let collectionMomentJson = []
        if (
          element[1].collectionMoment !== '' &&
          element[1].collectionMoment !== undefined
        ) {
          if (typeof element[1].collectionMoment === 'string') {
            collectionMomentJson = element[1].collectionMoment
          }
        }

        let obtentionProcedureJson = []

        if (
          element[1].obtentionProcedure !== '' &&
          element[1].obtentionProcedure !== undefined
        ) {
          if (typeof element[1].obtentionProcedure === 'object') {
            element[1].obtentionProcedure.forEach(element2 => {
              obtentionProcedureJson.push(
                JSON.stringify(element2, null, 2)
                  .replaceAll('[', '')
                  .replaceAll(']', '')
                  .replaceAll('{', '')
                  .replaceAll('}', '')
                  .replaceAll(',', '')
                  .replaceAll(' ,', '')
                  .replaceAll(', ', '')
                  .replaceAll('"', '')
              )
            })
            obtentionProcedureJson = obtentionProcedureJson.toString()
            obtentionProcedureJson = obtentionProcedureJson
              .replaceAll(', ', ',')
              .replaceAll(' ,', ',')
            obtentionProcedureJson = obtentionProcedureJson.replaceAll(',', '')
          } else {
            obtentionProcedureJson = JSON.stringify(
              element[1].obtentionProcedure,
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
            obtentionProcedureJson = obtentionProcedureJson.toString()
            obtentionProcedureJson = obtentionProcedureJson
              .replaceAll(', ', ',')
              .replaceAll(' ,', ',')
            obtentionProcedureJson = obtentionProcedureJson.replaceAll(',', '')
          }
        }

        let tumorProgressionJson = []

        if (
          element[1].tumorProgression !== '' &&
          element[1].tumorProgression !== undefined
        ) {
          if (typeof element[1].tumorProgression === 'object') {
            element[1].tumorProgression.forEach(element2 => {
              tumorProgressionJson.push(
                JSON.stringify(element2, null, 2)
                  .replaceAll('[', '')
                  .replaceAll(']', '')
                  .replaceAll('{', '')
                  .replaceAll('}', '')
                  .replaceAll(',', '')
                  .replaceAll(' ,', '')
                  .replaceAll(', ', '')
                  .replaceAll('"', '')
              )
            })
            tumorProgressionJson = tumorProgressionJson.toString()
            tumorProgressionJson = tumorProgressionJson
              .replaceAll(', ', ',')
              .replaceAll(' ,', ',')
            tumorProgressionJson = tumorProgressionJson.replaceAll(',', '')
          } else {
            tumorProgressionJson = JSON.stringify(
              element[1].tumorProgression,
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
            tumorProgressionJson = tumorProgressionJson.toString()
            tumorProgressionJson = tumorProgressionJson
              .replaceAll(', ', ',')
              .replaceAll(' ,', ',')
            tumorProgressionJson = tumorProgressionJson.replaceAll(',', '')
          }
        }

        let tumorGradeJson = []

        if (
          element[1].tumorGrade !== '' &&
          element[1].tumorGrade !== undefined
        ) {
          if (typeof element[1].tumorGrade === 'object') {
            element[1].tumorGrade.forEach(element2 => {
              tumorGradeJson.push(
                JSON.stringify(element2, null, 2)
                  .replaceAll('[', '')
                  .replaceAll(']', '')
                  .replaceAll('{', '')
                  .replaceAll('}', '')
                  .replaceAll(',', '')
                  .replaceAll(' ,', '')
                  .replaceAll(', ', '')
                  .replaceAll('"', '')
              )
            })
            tumorGradeJson = tumorGradeJson.toString()
            tumorGradeJson = tumorGradeJson
              .replaceAll(', ', ',')
              .replaceAll(' ,', ',')
            tumorGradeJson = tumorGradeJson.replaceAll(',', '')
          } else {
            tumorGradeJson = JSON.stringify(element[1].tumorGrade, null, 2)
              .replaceAll('[', '')
              .replaceAll(']', '')
              .replaceAll('{', '')
              .replaceAll('}', '')
              .replaceAll(',', '')
              .replaceAll(' ,', '')
              .replaceAll(', ', '')
              .replaceAll('"', '')
            tumorGradeJson = tumorGradeJson.toString()
            tumorGradeJson = tumorGradeJson
              .replaceAll(', ', ',')
              .replaceAll(' ,', ',')
            tumorGradeJson = tumorGradeJson.replaceAll(',', '')
          }
        }

        let pathologicalStageJson = {}

        if (
          element[1].pathologicalStage !== '' &&
          element[1].pathologicalStage !== undefined
        ) {
          if (typeof element[1].pathologicalStage === 'object') {
            pathologicalStageJson = {
              id: element[1].pathologicalStage.id,
              label: element[1].pathologicalStage.label
            }
          }
        }

        let pathologicalTnmFindingJson = []

        if (
          element[1].pathologicalTnmFinding !== '' &&
          element[1].pathologicalTnmFinding !== undefined
        ) {
          if (typeof element[1].pathologicalTnmFinding === 'object') {
            element[1].pathologicalTnmFinding.forEach(element2 => {
              pathologicalTnmFindingJson.push(
                JSON.stringify(element2, null, 2)
                  .replaceAll('[', '')
                  .replaceAll(']', '')
                  .replaceAll('{', '')
                  .replaceAll('}', '')
                  .replaceAll(',', '')
                  .replaceAll(' ,', '')
                  .replaceAll(', ', '')
                  .replaceAll('"', '')
              )
            })
            pathologicalTnmFindingJson = pathologicalTnmFindingJson.toString()
            pathologicalTnmFindingJson = pathologicalTnmFindingJson
              .replaceAll(', ', ',')
              .replaceAll(' ,', ',')
            pathologicalTnmFindingJson = pathologicalTnmFindingJson.replaceAll(
              ',',
              ''
            )
          } else {
            pathologicalTnmFindingJson = JSON.stringify(
              element[1].pathologicalTnmFinding,
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
            pathologicalTnmFindingJson = pathologicalTnmFindingJson.toString()
            pathologicalTnmFindingJson = pathologicalTnmFindingJson
              .replaceAll(', ', ',')
              .replaceAll(' ,', ',')
            pathologicalTnmFindingJson = pathologicalTnmFindingJson.replaceAll(
              ',',
              ''
            )
          }
        }

        let histologicalDiagnosis_id = ''
        let histologicalDiagnosis_label = ''
        let stringHistologicalDiagnosis = ''

        if (
          element[1].histologicalDiagnosis !== '' &&
          element[1].histologicalDiagnosis !== undefined
        ) {
          histologicalDiagnosis_id = element[1].histologicalDiagnosis.id
          histologicalDiagnosis_label = element[1].histologicalDiagnosis.label
          stringHistologicalDiagnosis = `${histologicalDiagnosis_id} / ${histologicalDiagnosis_label}`
        } else {
          stringHistologicalDiagnosis = ''
        }

        let diagnosticMarkersJson = []

        if (
          element[1].diagnosticMarkers !== '' &&
          element[1].diagnosticMarkers !== undefined
        ) {
          if (typeof element[1].diagnosticMarkers === 'object') {
            element[1].diagnosticMarkers.forEach(element2 => {
              diagnosticMarkersJson.push(
                JSON.stringify(element2, null, 2)
                  .replaceAll('[', '')
                  .replaceAll(']', '')
                  .replaceAll('{', '')
                  .replaceAll('}', '')
                  .replaceAll(',', '')
                  .replaceAll(' ,', '')
                  .replaceAll(', ', '')
                  .replaceAll('"', '')
              )
            })
            diagnosticMarkersJson = diagnosticMarkersJson.toString()
            diagnosticMarkersJson = diagnosticMarkersJson
              .replaceAll(', ', ',')
              .replaceAll(' ,', ',')
            diagnosticMarkersJson = diagnosticMarkersJson.replaceAll(',', '')
          } else {
            diagnosticMarkersJson = JSON.stringify(
              element[1].diagnosticMarkers,
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
            diagnosticMarkersJson = diagnosticMarkersJson.toString()
            diagnosticMarkersJson = diagnosticMarkersJson
              .replaceAll(', ', ',')
              .replaceAll(' ,', ',')
            diagnosticMarkersJson = diagnosticMarkersJson.replaceAll(',', '')
          }
        }

        let phenotypicFeaturesJson = []

        if (
          element[1].phenotypicFeatures !== '' &&
          element[1].phenotypicFeatures !== undefined
        ) {
          if (typeof element[1].phenotypicFeatures === 'object') {
            element[1].phenotypicFeatures.forEach(element2 => {
              phenotypicFeaturesJson.push(
                JSON.stringify(element2, null, 2)
                  .replaceAll('[', '')
                  .replaceAll(']', '')
                  .replaceAll('{', '')
                  .replaceAll('}', '')
                  .replaceAll(',', '')
                  .replaceAll(' ,', '')
                  .replaceAll(', ', '')
                  .replaceAll('"', '')
              )
            })
            phenotypicFeaturesJson = phenotypicFeaturesJson.toString()
            phenotypicFeaturesJson = phenotypicFeaturesJson
              .replaceAll(', ', ',')
              .replaceAll(' ,', ',')
            phenotypicFeaturesJson = phenotypicFeaturesJson.replaceAll(',', '')
          } else {
            phenotypicFeaturesJson = JSON.stringify(
              element[1].phenotypicFeatures,
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
            phenotypicFeaturesJson = phenotypicFeaturesJson.toString()
            phenotypicFeaturesJson = phenotypicFeaturesJson
              .replaceAll(', ', ',')
              .replaceAll(' ,', ',')
            phenotypicFeaturesJson = phenotypicFeaturesJson.replaceAll(',', '')
          }
        }

        let measurementsJson = []

        if (
          element[1].measurements !== '' &&
          element[1].measurements !== undefined
        ) {
          if (typeof element[1].measurements === 'object') {
            element[1].measurements.forEach(element2 => {
              measurementsJson.push(
                JSON.stringify(element2, null, 2)
                  .replaceAll('[', '')
                  .replaceAll(']', '')
                  .replaceAll('{', '')
                  .replaceAll('}', '')
                  .replaceAll(',', '')
                  .replaceAll(' ,', '')
                  .replaceAll(', ', '')
                  .replaceAll('"', '')
              )
            })
            measurementsJson = measurementsJson.toString()
            measurementsJson = measurementsJson
              .replaceAll(', ', ',')
              .replaceAll(' ,', ',')
            measurementsJson = measurementsJson.replaceAll(',', '')
          } else {
            measurementsJson = JSON.stringify(element[1].measurements, null, 2)
              .replaceAll('[', '')
              .replaceAll(']', '')
              .replaceAll('{', '')
              .replaceAll('}', '')
              .replaceAll(',', '')
              .replaceAll(' ,', '')
              .replaceAll(', ', '')
              .replaceAll('"', '')
            measurementsJson = measurementsJson.toString()
            measurementsJson = measurementsJson
              .replaceAll(', ', ',')
              .replaceAll(' ,', ',')
            measurementsJson = measurementsJson.replaceAll(',', '')
          }
        }

        let sampleProcessingJson = []

        if (
          element[1].sampleProcessing !== '' &&
          element[1].sampleProcessing !== undefined
        ) {
          if (typeof element[1].sampleProcessing === 'object') {
            element[1].sampleProcessing.forEach(element2 => {
              sampleProcessingJson.push(
                JSON.stringify(element2, null, 2)
                  .replaceAll('[', '')
                  .replaceAll(']', '')
                  .replaceAll('{', '')
                  .replaceAll('}', '')
                  .replaceAll(',', '')
                  .replaceAll(' ,', '')
                  .replaceAll(', ', '')
                  .replaceAll('"', '')
              )
            })
            sampleProcessingJson = sampleProcessingJson.toString()
            sampleProcessingJson = sampleProcessingJson
              .replaceAll(', ', ',')
              .replaceAll(' ,', ',')
            sampleProcessingJson = sampleProcessingJson.replaceAll(',', '')
          } else {
            sampleProcessingJson = JSON.stringify(
              element[1].sampleProcessing,
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
            sampleProcessingJson = sampleProcessingJson.toString()
            sampleProcessingJson = sampleProcessingJson
              .replaceAll(', ', ',')
              .replaceAll(' ,', ',')
            sampleProcessingJson = sampleProcessingJson.replaceAll(',', '')
          }
        }

        let sampleStorage_id = ''
        let sampleStorage_label = ''
        let stringSampleStorage = ''

        if (
          element[1].sampleStorage !== '' &&
          element[1].sampleStorage !== undefined
        ) {
          sampleStorage_id = element[1].sampleStorage.id
          sampleStorage_label = element[1].sampleStorage.label
          stringSampleStorage = `${sampleStorage_id} / ${sampleStorage_label}`
        } else {
          stringSampleStorage = ''
        }

        var myObjRows = new Object()
        myObjRows.id = index
        if (element[1].id !== '') {
          myObjRows.BiosampleId = element[1].id
        }
        if (element[1].individualId !== '') {
          myObjRows.individualId = element[1].individualId
        }
        myObjRows.Beacon = element[0]

        if (stringBiosampleStatus !== '') {
          myObjRows.biosampleStatus = stringBiosampleStatus
        }
        if (stringSampleOriginType !== '') {
          myObjRows.sampleOriginType = stringSampleOriginType
        }
        if (stringSampleOriginDetail !== '') {
          myObjRows.sampleOriginDetail = stringSampleOriginDetail
        }
        if (collectionDateJson !== '') {
          myObjRows.collectionDate = collectionDateJson
        }
        if (collectionMomentJson !== '') {
          myObjRows.collectionMoment = collectionMomentJson
        }
        if (obtentionProcedureJson !== '') {
          myObjRows.obtentionProcedure = obtentionProcedureJson
        }
        if (tumorProgressionJson !== '') {
          myObjRows.tumorProgression = tumorProgressionJson
        }
        if (tumorGradeJson !== '') {
          myObjRows.tumorGrade = tumorGradeJson
        }

        if (pathologicalStageJson !== '') {
          myObjRows.pathologicalStage = pathologicalStageJson
        }

        if (pathologicalTnmFindingJson !== '') {
          myObjRows.pathologicalTnmFinding = pathologicalTnmFindingJson
        }

        if (stringHistologicalDiagnosis !== '') {
          myObjRows.histologicalDiagnosis = stringHistologicalDiagnosis
        }
        if (diagnosticMarkersJson !== '') {
          myObjRows.diagnosticMarkers = diagnosticMarkersJson
        }
        if (phenotypicFeaturesJson !== '') {
          myObjRows.phenotypicFeatures = phenotypicFeaturesJson
        }
        if (measurementsJson !== '') {
          myObjRows.measurements = measurementsJson
        }
        if (sampleProcessingJson !== '') {
          myObjRows.sampleProcessing = sampleProcessingJson
        }
        if (stringSampleStorage !== '') {
          myObjRows.sampleStorage = stringSampleStorage
        }
        console.log(rows)
        rows.push(myObjRows)

        if (index === resultsSelectedFinal.length - 1) {
          setEditable(rows.map(o => ({ ...o })))

          setTrigger2(true)
        }
      }
    })
  }, [trigger, resultsSelectedFinal])

  useEffect(() => {
    let count = 0

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
          <DataGrid
            getRowHeight={() => 'auto'}
            checkboxSelection
            columns={columns}
            rows={editable}
            slots={{ toolbar: CustomToolbar }}
            slotProps={{
              toolbar: {
                printOptions: { getRowsToExport: getSelectedRowsToExport }
              }
            }}
          />
        )}
      {showCrossQuery && (
        <CrossQueries
          parameter={parameterCrossQuery}
          collection={'biosamples'}
          setShowCrossQuery={setShowCrossQuery}
        />
      )}
    </div>
  )
}

export default TableResultsBiosamples
