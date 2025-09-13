# PPI Affinity Predictor

A machine learning-based tool for predicting protein-protein interaction (PPI) binding affinity using protein sequences.

## Features

- **Sequence-based prediction**: Input two protein sequences and get predicted binding affinity
- **Multiple input formats**: Support for direct sequence input or FASTA files
- **Comprehensive feature extraction**: 
  - Amino acid composition
  - Dipeptide composition
  - Physicochemical properties
  - Structural features
- **Interactive mode**: User-friendly command-line interface
- **Model persistence**: Save and load trained models

## Installation

### Prerequisites

- Python 3.7 or higher
- pip package manager

### Setup Environment

1. Clone or download this repository
2. Navigate to the PPI_detect directory
3. Install dependencies:

```bash
pip install -r requirements.txt
```

### Alternative: Using conda

```bash
conda create -n ppi_env python=3.8
conda activate ppi_env
pip install -r requirements.txt
```

## Usage

### Quick Start

1. **Train a model** (first time usage):
```bash
python ppi_affinity_predictor.py --train
```

2. **Interactive mode** (recommended for beginners):
```bash
python ppi_affinity_predictor.py --interactive
```

3. **Command line prediction**:
```bash
python ppi_affinity_predictor.py --seq1 "MKFLILLFNILCLFPVLAADNHGVGPQGASLFRLAKEGCMVVLLLLLLLLLLLLLLLLLLLLFRLVKRLRAKEGCMVVLLLLLLLLLLLLLL" --seq2 "MSTLPLLLLLLLLPLHRAADNHGVGPQGASLFRLAKEGCMV"
```

### Detailed Usage Options

#### Training a New Model
```bash
python ppi_affinity_predictor.py --train --model my_model.pkl
```

#### Using FASTA Files
```bash
python ppi_affinity_predictor.py --seq1 protein1.fasta --seq2 protein2.fasta
```

#### Interactive Mode
```bash
python ppi_affinity_predictor.py --interactive --model my_model.pkl
```

### Command Line Arguments

- `--seq1`: First protein sequence (string or FASTA file path)
- `--seq2`: Second protein sequence (string or FASTA file path)
- `--train`: Train a new model with synthetic data
- `--model`: Path to model file (default: ppi_model.pkl)
- `--interactive`: Run in interactive mode

## Input Format

### Protein Sequences
- **Direct input**: Single-letter amino acid codes (e.g., "MKFLILLFNIL...")
- **FASTA files**: Standard FASTA format files

### Example FASTA File
```
>Protein1
MKFLILLFNILCLFPVLAADNHGVGPQGASLFRLAKEGCMVVLLLLLLLLLLLLL
LLLLLLFRLVKRLRAKEGCMVVLLLLLLLLLLLLLL
>Protein2
MSTLPLLLLLLLLPLHRAADNHGVGPQGASLFRLAKEGCMV
```

## Output

The predictor outputs binding affinity values in pKd units:

- **pKd > 9**: Very strong binding
- **pKd 7-9**: Strong binding  
- **pKd 5-7**: Moderate binding
- **pKd < 5**: Weak binding

### Example Output
```
Predicted binding affinity: 7.245 (pKd)
Binding strength interpretation:
- Strong binding
```

## Features Extracted

The predictor analyzes various protein sequence features:

### Compositional Features
- Amino acid composition (20 features)
- Dipeptide composition (400 features)

### Physicochemical Properties
- Molecular weight
- Isoelectric point
- Instability index
- Flexibility
- Aromaticity
- GRAVY (hydropathy)
- Secondary structure fractions

### Structural Features
- Sequence length
- Hydrophobic residue ratio
- Polar residue ratio
- Charged residue ratio
- Aromatic residue ratio

### Interaction Features
- Length ratios and differences
- Feature differences between proteins
- Feature ratios between proteins

## Model Architecture

- **Algorithm**: Random Forest Regression
- **Features**: ~450+ sequence-derived features
- **Preprocessing**: StandardScaler normalization
- **Validation**: Train-test split with performance metrics

## Performance Metrics

The model reports:
- **MSE**: Mean Squared Error
- **R²**: Coefficient of determination

## Example Scripts

### Simple Prediction
```python
from ppi_affinity_predictor import PPIAffinityPredictor

predictor = PPIAffinityPredictor()
predictor.load_model('ppi_model.pkl')

seq1 = "MKFLILLFNILCLFPVLA"
seq2 = "MSTLPLLLLLLLLPLHRA" 

affinity = predictor.predict_affinity(seq1, seq2)
print(f"Predicted affinity: {affinity:.3f}")
```

### Batch Prediction
```python
sequences_pairs = [
    ("MKFLILLFNIL...", "MSTLPLLLLLL..."),
    ("ATGCGTACGTA...", "CCGTAGCTAAG..."),
]

for seq1, seq2 in sequences_pairs:
    affinity = predictor.predict_affinity(seq1, seq2)
    print(f"Affinity: {affinity:.3f}")
```

## Troubleshooting

### Common Issues

1. **ModuleNotFoundError**: Install missing dependencies
   ```bash
   pip install -r requirements.txt
   ```

2. **Invalid sequence characters**: Ensure sequences contain only standard amino acids
   - Valid: A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y

3. **Memory issues**: For very long sequences (>1000 residues), consider splitting or using a more powerful machine

4. **Model not found**: Train a new model first:
   ```bash
   python ppi_affinity_predictor.py --train
   ```

### Performance Tips

- Use shorter sequences when possible (50-500 residues work best)
- For batch processing, load the model once and reuse
- Consider using GPU acceleration for larger datasets

## Limitations

- **Synthetic training data**: Current model uses synthetic data for demonstration
- **Sequence length**: Optimized for sequences 50-500 residues
- **No structural information**: Based purely on sequence features
- **Species-specific**: May need retraining for specific organisms

## Future Improvements

- Integration with real experimental PPI affinity data
- Deep learning models (CNN, RNN, Transformers)
- Structural feature integration
- Species-specific models
- Web interface development

## Citation

If you use this tool in your research, please cite:

```
PPI Affinity Predictor
Author: YuGe0115
Year: 2025
GitHub: [Repository URL]
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Support

For questions or issues:
1. Check the troubleshooting section
2. Review the example scripts
3. Open an issue on GitHub

## Changelog

### Version 1.0.0 (2025)
- Initial release
- Basic sequence-based affinity prediction
- Interactive and command-line modes
- Synthetic data training pipeline