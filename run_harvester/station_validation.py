"""Station-set validation and collision-resistant log naming for ADCIRC jobs."""

import re

import numpy as np
import pandas as pd


def _station_ids(values):
    return {str(value).strip() for value in values if str(value).strip()}


def requested_station_ids(station_list):
    """Return station IDs from either a Fort61 list or Fort63 DataFrame."""
    if isinstance(station_list, pd.DataFrame):
        if 'stationid' not in station_list.columns:
            return set()
        return _station_ids(station_list['stationid'].tolist())
    return _station_ids(station_list)


def validate_station_sets(data, metadata, requested=None, policy='fail'):
    """Compare requested, metadata, and emitted data stations before CSV output.

    Metadata is the source-confirmed station set. Requested stations absent from
    both source metadata and data are reported but are not treated as an output
    mismatch because Fort61 station lists may intentionally be supersets.
    Source-confirmed all-fill columns are reported as data warnings because dry
    ADCIRC nodes are valid, but they are not structural station-set mismatches.
    """
    if policy not in {'fail', 'warn'}:
        raise ValueError(f'Unknown station mismatch policy: {policy}')

    requested_ids=_station_ids(requested or [])
    data_ids=_station_ids(data.columns)
    metadata_ids=_station_ids(metadata.index)
    no_valid_data=[]
    for station in sorted(data_ids):
        values=pd.to_numeric(data[station], errors='coerce')
        valid=np.isfinite(values.to_numpy(dtype=float)) & (values.to_numpy(dtype=float) != -99999)
        if not bool(valid.any()):
            no_valid_data.append(station)

    report={
        'requested': sorted(requested_ids),
        'metadata': sorted(metadata_ids),
        'data': sorted(data_ids),
        'requested_without_metadata': sorted(requested_ids-metadata_ids),
        'metadata_without_data': sorted(metadata_ids-data_ids),
        'data_without_metadata': sorted(data_ids-metadata_ids),
        'stations_without_valid_data': no_valid_data,
    }
    # A column that is present in both data and metadata can legitimately be
    # all fill (for example, an ADCIRC station on a dry node). Keep that
    # condition visible without treating it as structural station loss.
    report['has_mismatch']=bool(
        report['metadata_without_data']
        or report['data_without_metadata']
    )
    report['has_data_warning']=bool(report['stations_without_valid_data'])
    if report['has_mismatch'] and policy == 'fail':
        raise RuntimeError(f'ADCIRC station-set validation failed: {report}')
    return report


def safe_log_token(value, default):
    token=str(value or default).strip()
    return re.sub(r'[^A-Za-z0-9_.-]+', '-', token) or default


def log_identity(model_id, provider, cast, pod):
    return '__'.join(
        [
            safe_log_token(model_id, 'nomodel'),
            safe_log_token(provider, 'noprovider'),
            safe_log_token(cast, 'nocast'),
            safe_log_token(pod, 'nopod'),
        ]
    )


def archived_log_name(model_id, provider, cast, pod):
    return f'run_harvester_adcirc_{log_identity(model_id, provider, cast, pod)}.log'
