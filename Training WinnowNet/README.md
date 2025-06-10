# Training WinnowNet Models

This folder contains scripts, datasets, and instructions for training two variants of the WinnowNet deep learning model: a self-attention-based model and a CNN-based model. Training is carried out in two phases to enable curriculum learning from synthetic (easy) to real-world metaproteomic (difficult) datasets.

## Requirements

- Python 3.7+
- PyTorch
- NumPy, Pandas, scikit-learn

## Datasets

- **Prosit_train.zip** (Phase 1 training set):   https://figshare.com/articles/dataset/Datasets/25511770?file=55257041
- **marine1_train.zip** (Phase 2 training set): https://figshare.com/articles/dataset/Datasets/25511770?file=55257035

---

## Self-Attention-Based WinnowNet

### Phase 1: Training on Easy Tasks (Synthetic Data)

```bash
python SpectraFeatures_training.py -i filename.tsv -s filename.ms2 -o spectra_feature.pkl -t 20 -f att
python WinnowNet_Att.py -i spectra_feature_directory -m prosit_att.pt
```

**Explanation of options:**
- `-i`: Input tab-delimited file with PSMs, including labels and weights.
- `-s`: Corresponding MS2 file (filename should match TSV).
- `-o`: Output file to store extracted features as a `pkl` file.
- `-t`: Number of threads for parallel processing.
- `-f`: Feature type (`att` for self-attention model).
- `-m`: Filename to save the trained model.
- A for-loop is needed to convert all `tsv` files to `pkl` files.

### Phase 2: Training on Difficult Tasks (Real Data)

```bash
python SpectraFeatures_training.py -i filename.tsv -s filename.ms2 -o spectra_feature.pkl -t 20 -f att
python WinnowNet_Att.py -i spectra_feature_directory -m marine_att.pt -p prosit_att.pt
```

- `-p`: Pre-trained model from Phase 1.
- A for-loop is needed to convert all `tsv` files to `pkl` files.

**Pre-trained model:** marine_att.pt,  https://figshare.com/articles/dataset/Models/25513531

---

## CNN-Based WinnowNet

### Phase 1: Training on Easy Tasks (Synthetic Data)

```bash
python SpectraFeatures_training.py -i filename.tsv -s filename.ms2 -o spectra_feature.pkl -t 20 -f cnn
python WinnowNet_CNN.py -i spectra_feature_directory -m prosit_cnn.pt
```

### Phase 2: Training on Difficult Tasks (Real Data)

```bash
python SpectraFeatures_training.py -i filename.tsv -s filename.ms2 -o spectra_feature.pkl -t 20 -f cnn
python WinnowNet_Att.py -i spectra_feature_directory -m cnn_pytorch.pt -p prosit_cnn.pt
```

**Pre-trained model:** cnn_pytorch.pt, https://figshare.com/articles/dataset/Models/25513531

---

## Notes

- All input MS2/TSV files must be preprocessed properly.
- Models trained in Phase 1 are reused to initialize weights in Phase 2.
- Training with GPU is recommended for performance.
