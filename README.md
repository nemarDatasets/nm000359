# Bern-Barcelona EEG database (iEEG-BIDS)

Focal and non-focal intracranial EEG signal pairs from five patients with pharmacoresistant temporal lobe
epilepsy, published with:

> Andrzejak RG, Schindler K, Rummel C (2012). Nonrandomness, nonlinear dependence, and nonstationarity of
> electroencephalographic recordings from epilepsy patients. *Phys. Rev. E* 86, 046206.
> doi:[10.1103/PhysRevE.86.046206](https://doi.org/10.1103/PhysRevE.86.046206)

Source dataset: Andrzejak RG, Schindler KA, Rummel C. *Nonrandomness, nonlinear dependence, and nonstationarity
of electroencephalographic recordings from epilepsy patients [dataset]*. Repositori Digital de la UPF, 2012.
doi:[10.34810/data502](https://doi.org/10.34810/data502), hdl:[10230/42829](http://hdl.handle.net/10230/42829).
The authors ask users to refer to the resources as the "Bern-Barcelona EEG database" and to cite the paper.

This BIDS dataset is a lossless re-packaging of that source. Every sample is identical to the published text
files; nothing was filtered, resampled, re-referenced or rescaled during conversion.

## Content

| | |
|---|---|
| Signal pairs | 7500 (3750 focal + 3750 non-focal) |
| Channels per run | 2 (`x`, `y`: neighbouring channels of the same class) |
| Samples per channel | 10240 |
| Sampling rate | 512 Hz |
| Duration per run | 20 s |
| Patients | 5 (pooled, not identifiable per run; see below) |
| Value range | -2223.068115 to 3469.565918 (source units, see "Units") |

## Recordings (from the paper, Sec. II)

- Five patients with longstanding pharmacoresistant temporal lobe epilepsy, candidates for epilepsy surgery,
  underwent long-term intracranial EEG recordings at the Department of Neurology of the University of Bern
  (Inselspital) because non-invasive studies had not localized the seizure onset zone.
- Intracranial strip and depth electrodes, all manufactured by AD-TECH (Racine, WI, USA). An extracranial
  reference electrode was placed between 10-20 positions Fz and Pz.
- Signals were sampled at 512 or 1024 Hz depending on whether they were recorded with more or less than 64
  channels.
- All five patients had good surgical outcome (three seizure-free; two with auras only; ILAE class 1 and 2).
- Ethics: "Retrospective EEG data analysis has been approved by the ethics committee of the Kanton of Bern. In
  addition, all patients gave written informed consent that their data from long-term EEG might be used for
  research purposes."

## Preprocessing applied by the source (before publication)

From the paper, Sec. II B, and the source code comments (`ASR_Setparameters.m`):

1. Digital band-pass filter 0.5-150 Hz, fourth-order Butterworth, applied forward and backward (zero phase).
2. Signals recorded at 1024 Hz were down-sampled to 512 Hz. Which pairs were down-sampled is not recorded.
3. Re-referencing against the median of all channels free of permanent artifacts (visual inspection).

`ASR_Setparameters.m`: "The data as provided on the page is sampled at 512 HZ and band-pass filtered between
0.5Hz and 150Hz." The additional 40 Hz low-pass and 46.5-53.5 Hz band-stop filters used in the paper's analysis
are part of the analysis code. They were not applied to the published data, and they were not applied here.

## Selection of the signal pairs (paper, Sec. II C)

- "Focal EEG channels": all channels that detected the first ictal EEG signal changes, as judged by visual
  inspection by at least two board-certified electroencephalographers. All other channels are "non-focal".
- The recordings were divided into 20-s windows (10240 samples). Seizure recordings and the three hours after
  the last seizure were excluded.
- For each focal pair the authors randomly selected a patient, one of that patient's focal channels (signal
  `x`), one neighbouring focal channel (signal `y`) and one time window. Sampling was uniform and without
  replacement. Non-focal pairs were drawn the same way from non-focal channels.
- Pairs with prominent measurement artifacts were discarded after visual inspection. Moderate 50 Hz line noise
  was not an exclusion criterion. No clinical criteria (presence or absence of epileptiform activity) were
  applied.
- "the focal EEG signal pairs were stored in the order in which they were drawn. Their origin (patient,
  channel, window) was not stored."

## BIDS layout and source-to-BIDS mapping

The source holds no patient, channel or time information for any pair, so per-patient subjects cannot be
built. All runs sit under one pseudo-subject, `sub-pooled`. This is **not one person**: it pools anonymous
segments from five patients. Do not treat runs as independent subjects. Pairs from the same patient cannot
be grouped either, so cross-validation by patient is not possible with this dataset.

| Source | BIDS |
|---|---|
| `Data_F_Ind_<a>_<b>.zip` / `Data_F_Ind<NNNN>.txt` (focal pair NNNN) | `sub-pooled/ieeg/sub-pooled_task-interictal_acq-focal_run-<NNNN>_ieeg.{vhdr,vmrk,eeg}` |
| `Data_N_Ind_<a>_<b>.zip` / `Data_N_Ind<NNNN>.txt` (non-focal pair NNNN) | `sub-pooled/ieeg/sub-pooled_task-interictal_acq-nonfocal_run-<NNNN>_ieeg.{vhdr,vmrk,eeg}` |
| column 1 / column 2 of each text file | channels `x` / `y` |
| `Results_F_All.txt`, `Results_N_All.txt` (published test outcomes) | columns of `sub-pooled/sub-pooled_scans.tsv` (described in `sub-pooled_scans.json`) |
| `Data_F_50.zip`, `Data_N_50.zip` (preview subsets) | not converted: their 100 members are byte-identical (SHA-256) to pairs 1-50 of the full archives |
| `ASR_Sources_2013_06_11.zip` (MATLAB analysis code) | kept in `sourcedata/` only |
| all of the above + repository metadata | byte-identical copies in `sourcedata/upf-repositori-10230-42829/` |

The run number equals the source pair index (`run-0007` = `Data_F_Ind0007.txt` for `acq-focal`).
`sub-pooled_scans.tsv` lists for every run the source archive, source file, the SHA-256 of the
uncompressed source file and the published test outcomes. Sidecars are shared by inheritance:
`sub-pooled_task-interictal_ieeg.json`, `sub-pooled_task-interictal_channels.tsv`, and one
`..._acq-<focal|nonfocal>_events.tsv` per class with a single event spanning the run (`trial_type` =
`focal` or `nonfocal`).

`task-interictal` is a label for the brain state (seizure-free interval of clinical monitoring). There
was no task. `acq-focal` / `acq-nonfocal` encode the signal class defined by the source.

## Data format and exactness

The source stores values as text with six decimals. Every value was converted to IEEE float32 (BrainVision
`IEEE_FLOAT_32`, multiplexed, resolution 1). For all 7500 files and all 153,600,000 values, formatting the
stored float32 value with `%.6f` reproduces the source text token exactly. The conversion is therefore
lossless with respect to the published text. A separate round-trip check re-read every BrainVision file and
compared it with the source text (see the campaign ledger).

## Units

Neither the paper nor the download page states the physical unit or scale of the stored values. The paper's
example figure shows values of roughly +/-100 to +/-200 with no unit label. Channel `units` are therefore
`n/a` in `channels.tsv` and in the BrainVision header. Values are exactly the source numbers. Magnitudes are
consistent with microvolts, but this is not documented by the source.

Channel `type` is `OTHER` because the source does not record whether a signal comes from a strip (ECoG) or a
depth (SEEG) contact. Electrode positions are not available, so there is no `electrodes.tsv`.

## Privacy

The published files contain only numeric samples (no headers, names, dates or identifiers). The source
randomized pairs across patients and did not store patient, channel or time of origin. Converted headers
contain no dates or identifiers. Repository metadata (`item.json`, `bundles.json`) contains only
bibliographic information.

## License and terms of use

- Repository record (dc.rights): "Licensed under a Creative Commons License (CC-BY) 4.0"
  (https://creativecommons.org/licenses/by/4.0/). This dataset is released under **CC-BY-4.0**.
- The UPF NTSA download page for the same material (https://www.upf.edu/web/ntsa/downloads, captured
  2026-10-06; text in `sourcedata/upf-repositori-10230-42829/_repository_metadata/`) adds a "Legal Agreement":
  "The source codes, data and results on these sites are free of charge for research and education purposes
  only. Any commercial or military use is prohibited. All resources are provided without any expressed or
  implied warranty. In no event the authors of the article or any of their host institutions are liable for any
  damages arising from the use of the software, data or results."
  This statement is not part of the repository record's rights field for this item. It is reproduced here so
  users can respect the authors' stated intent.

## Funding

R.G.A. acknowledges grant FIS-2010-18204 of the Spanish Ministry of Education and Science.

## How to load

Example with MNE-Python / MNE-BIDS (after `nemar dataset get` or `datalad get` of the files you need):

```python
import mne
from mne_bids import BIDSPath, read_raw_bids

bp = BIDSPath(root="<dataset root>", subject="pooled", task="interictal",
              acquisition="focal", run="0007", datatype="ieeg", suffix="ieeg", extension=".vhdr")
raw = read_raw_bids(bp)          # 2 channels (x, y), 10240 samples at 512 Hz
data = raw.get_data()            # values exactly as in the source text; physical unit not stated by the source
```

The published test outcomes per pair are in
`sub-pooled/sub-pooled_scans.tsv` (e.g. `pandas.read_csv(..., sep="\t")`).

## How to cite

Cite the paper (doi:10.1103/PhysRevE.86.046206) and the dataset (doi:10.34810/data502). Also cite this BIDS
release by its NEMAR identifier.

## References

- Andrzejak RG, Schindler K, Rummel C (2012). Phys. Rev. E 86, 046206. doi:10.1103/PhysRevE.86.046206
- Andrzejak RG, Schindler KA, Rummel C (2012). Dataset, Repositori Digital de la UPF. doi:10.34810/data502
- Andrzejak RG et al. (2001). Phys. Rev. E 64, 061907. doi:10.1103/PhysRevE.64.061907 (methods referenced by
  the analysis code)

## Conversion provenance

Converted 2026-10-06 on SDSC Voyager (Kubernetes jobs) by the iEEG-NEMAR campaign (lane G) with
`laneG_convert.py` and `laneG_finalize.py`. Source files were downloaded through the repository's DSpace REST
API, and every bitstream's MD5 matched the repository checksum. See `sourcedata/provenance.json`. Paper
details were taken from the published version deposited at hdl:10230/43557 (repository full-text extraction)
and from the UPF NTSA download page.
