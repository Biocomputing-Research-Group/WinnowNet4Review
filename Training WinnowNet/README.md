# Training WinnowNet Models

This folder contains the scripts and instructions for training the WinnowNet models used in our research. Two architectures are supported: a self-attention-based model and a CNN-based model. The training follows a two-phase strategy to ensure effective curriculum learning from synthetic to real-world metaproteomic datasets.

## Requirements

Before running the scripts, ensure the following dependencies are installed:
- Python 3.7+
- PyTorch
- NumPy, Pandas, scikit-learn
- Pickle (for handling `.pkl` feature files)

## Data Sources

- **Prosit training PSMs** (Phase 1): [Prosit_training_PSMs.tsv](https://figshare.com/articles/dataset/Datasets/25511770?file=54794327)
- **Marine1 training PSMs** (Phase 2): [marine1_training.tsv](https://figshare.com/articles/dataset/Datasets/25511770?file=54794327)

---

## Self-Attention-Based WinnowNet

### Phase 1 Training

```
python WinnowNet_Att.py -i spectra_feature.pkl -m prosit_att.pt
```

- `-i spectra_feature.pkl`: Input feature file extracted for each PSM
- `-m prosit_att.pt`: Output model filename for saving the trained model

### Phase 2 Training

```
python WinnowNet_Att.py -i spectra_feature.pkl -m marine_att.pt -p prosit_att.pt
```

- `-p prosit_att.pt`: Pre-trained model from Phase 1 used to initialize weights

**Pre-trained model download:** [marine_att.pt](https://figshare.com/articles/dataset/Models/25513531)

---

## CNN-Based WinnowNet

### Phase 1 Training

```
python WinnowNet_CNN.py -i spectra_feature.pkl -m prosit_cnn.pt
```

- Similar arguments as above; this trains the CNN-based version.

### Phase 2 Training

```
python WinnowNet_Att.py -i spectra_feature.pkl -m cnn_pytorch.pt -p prosit_cnn.pt
```

- Reuses the same script structure with the CNN model for transfer learning.

**Pre-trained model download:** [cnn_pytorch.pt](https://figshare.com/articles/dataset/Models/25513531)

---

## Notes

- Ensure all `.pkl` feature files are preprocessed prior to training.
- Training can be done using GPU for faster performance.

