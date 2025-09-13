#!/usr/bin/env python3
"""
Protein-Protein Interaction Affinity Predictor

This script predicts the binding affinity between two protein sequences
using various features including amino acid composition, physicochemical properties,
and sequence-based descriptors.
"""

import numpy as np
import pandas as pd
import argparse
import sys
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
import pickle
import os
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from itertools import product
import warnings
warnings.filterwarnings('ignore')


class ProteinFeatureExtractor:
    """Extract various features from protein sequences"""
    
    def __init__(self):
        self.amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
        self.hydrophobic = 'AILMFPWVY'
        self.polar = 'NQST'
        self.charged = 'DEKR'
        self.aromatic = 'FWY'
        
    def amino_acid_composition(self, sequence):
        """Calculate amino acid composition"""
        sequence = sequence.upper()
        composition = {}
        for aa in self.amino_acids:
            composition[f'AA_{aa}'] = sequence.count(aa) / len(sequence)
        return composition
    
    def dipeptide_composition(self, sequence):
        """Calculate dipeptide composition"""
        sequence = sequence.upper()
        dipeptides = [''.join(p) for p in product(self.amino_acids, repeat=2)]
        composition = {}
        
        if len(sequence) < 2:
            for dp in dipeptides:
                composition[f'DP_{dp}'] = 0.0
            return composition
            
        total_dipeptides = len(sequence) - 1
        for dp in dipeptides:
            count = sum(1 for i in range(len(sequence) - 1) 
                       if sequence[i:i+2] == dp)
            composition[f'DP_{dp}'] = count / total_dipeptides
        return composition
    
    def physicochemical_properties(self, sequence):
        """Calculate physicochemical properties"""
        try:
            analyzer = ProteinAnalysis(sequence.upper())
            properties = {}
            
            properties['molecular_weight'] = analyzer.molecular_weight()
            properties['isoelectric_point'] = analyzer.isoelectric_point()
            properties['instability_index'] = analyzer.instability_index()
            properties['flexibility'] = np.mean(analyzer.flexibility())
            properties['aromaticity'] = analyzer.aromaticity()
            properties['gravy'] = analyzer.gravy()
            
            sec_struct = analyzer.secondary_structure_fraction()
            properties['helix_fraction'] = sec_struct[0]
            properties['turn_fraction'] = sec_struct[1]
            properties['sheet_fraction'] = sec_struct[2]
            
        except Exception as e:
            print(f"Warning: Error calculating properties for sequence: {e}")
            properties = {
                'molecular_weight': 0, 'isoelectric_point': 7,
                'instability_index': 0, 'flexibility': 0,
                'aromaticity': 0, 'gravy': 0,
                'helix_fraction': 0, 'turn_fraction': 0, 'sheet_fraction': 0
            }
        
        return properties
    
    def structural_features(self, sequence):
        """Calculate structural features"""
        sequence = sequence.upper()
        length = len(sequence)
        
        features = {
            'length': length,
            'hydrophobic_ratio': sum(1 for aa in sequence if aa in self.hydrophobic) / length,
            'polar_ratio': sum(1 for aa in sequence if aa in self.polar) / length,
            'charged_ratio': sum(1 for aa in sequence if aa in self.charged) / length,
            'aromatic_ratio': sum(1 for aa in sequence if aa in self.aromatic) / length,
        }
        
        return features
    
    def extract_features(self, sequence):
        """Extract all features from a protein sequence"""
        features = {}
        features.update(self.amino_acid_composition(sequence))
        features.update(self.physicochemical_properties(sequence))
        features.update(self.structural_features(sequence))
        
        return features


class PPIAffinityPredictor:
    """Protein-Protein Interaction Affinity Predictor"""
    
    def __init__(self):
        self.feature_extractor = ProteinFeatureExtractor()
        self.model = None
        self.scaler = None
        self.feature_names = None
        
    def extract_ppi_features(self, seq1, seq2):
        """Extract features for protein-protein interaction"""
        features1 = self.feature_extractor.extract_features(seq1)
        features2 = self.feature_extractor.extract_features(seq2)
        
        ppi_features = {}
        
        # Individual protein features
        for key, value in features1.items():
            ppi_features[f'prot1_{key}'] = value
        for key, value in features2.items():
            ppi_features[f'prot2_{key}'] = value
            
        # Interaction features
        ppi_features['length_ratio'] = len(seq1) / len(seq2) if len(seq2) > 0 else 1
        ppi_features['length_diff'] = abs(len(seq1) - len(seq2))
        ppi_features['length_sum'] = len(seq1) + len(seq2)
        
        # Feature differences and ratios
        for key in features1.keys():
            if features1[key] != 0 and features2[key] != 0:
                ppi_features[f'ratio_{key}'] = features1[key] / features2[key]
            ppi_features[f'diff_{key}'] = abs(features1[key] - features2[key])
        
        return ppi_features
    
    def generate_synthetic_data(self, n_samples=1000):
        """Generate synthetic training data for demonstration"""
        print("Generating synthetic training data...")
        
        # Generate random protein sequences
        amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
        sequences1 = []
        sequences2 = []
        affinities = []
        
        np.random.seed(42)
        
        for _ in range(n_samples):
            # Generate sequences of varying lengths
            len1 = np.random.randint(50, 500)
            len2 = np.random.randint(50, 500)
            
            seq1 = ''.join(np.random.choice(list(amino_acids), len1))
            seq2 = ''.join(np.random.choice(list(amino_acids), len2))
            
            sequences1.append(seq1)
            sequences2.append(seq2)
            
            # Generate synthetic affinity based on simple rules
            # (In reality, this would be experimental data)
            length_factor = min(len1, len2) / 100
            hydrophobic_factor = (seq1.count('A') + seq1.count('L') + 
                                seq2.count('A') + seq2.count('L')) / (len1 + len2)
            charged_factor = (seq1.count('K') + seq1.count('R') + seq1.count('D') + 
                            seq2.count('K') + seq2.count('R') + seq2.count('D')) / (len1 + len2)
            
            # Synthetic affinity (pKd values typically range from 4-12)
            affinity = 6 + 2 * length_factor + 3 * hydrophobic_factor + 2 * charged_factor + np.random.normal(0, 0.5)
            affinity = max(4, min(12, affinity))  # Clamp to realistic range
            
            affinities.append(affinity)
        
        return sequences1, sequences2, affinities
    
    def train_model(self, sequences1, sequences2, affinities):
        """Train the affinity prediction model"""
        print("Extracting features from training data...")
        features_list = []
        
        for seq1, seq2 in zip(sequences1, sequences2):
            features = self.extract_ppi_features(seq1, seq2)
            features_list.append(features)
        
        # Convert to DataFrame
        df = pd.DataFrame(features_list)
        self.feature_names = df.columns.tolist()
        
        # Handle missing values
        df = df.fillna(0)
        
        X = df.values
        y = np.array(affinities)
        
        # Split data
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.2, random_state=42
        )
        
        # Scale features
        self.scaler = StandardScaler()
        X_train_scaled = self.scaler.fit_transform(X_train)
        X_test_scaled = self.scaler.transform(X_test)
        
        # Train model
        print("Training Random Forest model...")
        self.model = RandomForestRegressor(
            n_estimators=100,
            max_depth=10,
            random_state=42,
            n_jobs=-1
        )
        self.model.fit(X_train_scaled, y_train)
        
        # Evaluate model
        y_pred = self.model.predict(X_test_scaled)
        mse = mean_squared_error(y_test, y_pred)
        r2 = r2_score(y_test, y_pred)
        
        print(f"Model performance:")
        print(f"MSE: {mse:.3f}")
        print(f"R²: {r2:.3f}")
        
        return self.model
    
    def predict_affinity(self, seq1, seq2):
        """Predict binding affinity between two protein sequences"""
        if self.model is None:
            raise ValueError("Model not trained. Please train the model first.")
        
        # Extract features
        features = self.extract_ppi_features(seq1, seq2)
        
        # Convert to DataFrame and ensure feature order
        df = pd.DataFrame([features])
        df = df.reindex(columns=self.feature_names, fill_value=0)
        
        # Scale features
        X = self.scaler.transform(df.values)
        
        # Predict
        affinity = self.model.predict(X)[0]
        
        return affinity
    
    def save_model(self, filepath):
        """Save the trained model"""
        model_data = {
            'model': self.model,
            'scaler': self.scaler,
            'feature_names': self.feature_names
        }
        with open(filepath, 'wb') as f:
            pickle.dump(model_data, f)
        print(f"Model saved to {filepath}")
    
    def load_model(self, filepath):
        """Load a trained model"""
        with open(filepath, 'rb') as f:
            model_data = pickle.load(f)
        
        self.model = model_data['model']
        self.scaler = model_data['scaler']
        self.feature_names = model_data['feature_names']
        print(f"Model loaded from {filepath}")


def read_fasta(filepath):
    """Read sequences from FASTA file"""
    sequences = []
    for record in SeqIO.parse(filepath, "fasta"):
        sequences.append(str(record.seq))
    return sequences


def main():
    parser = argparse.ArgumentParser(
        description="Predict protein-protein interaction affinity"
    )
    parser.add_argument(
        "--seq1", type=str,
        help="First protein sequence (or path to FASTA file)"
    )
    parser.add_argument(
        "--seq2", type=str,
        help="Second protein sequence (or path to FASTA file)"
    )
    parser.add_argument(
        "--train", action="store_true",
        help="Train a new model with synthetic data"
    )
    parser.add_argument(
        "--model", type=str, default="ppi_model.pkl",
        help="Path to model file"
    )
    parser.add_argument(
        "--interactive", action="store_true",
        help="Run in interactive mode"
    )
    
    args = parser.parse_args()
    
    predictor = PPIAffinityPredictor()
    
    # Training mode
    if args.train:
        print("Training new model...")
        seq1_list, seq2_list, affinities = predictor.generate_synthetic_data()
        predictor.train_model(seq1_list, seq2_list, affinities)
        predictor.save_model(args.model)
        return
    
    # Load existing model
    if os.path.exists(args.model):
        predictor.load_model(args.model)
    else:
        print(f"Model file {args.model} not found. Training new model...")
        seq1_list, seq2_list, affinities = predictor.generate_synthetic_data()
        predictor.train_model(seq1_list, seq2_list, affinities)
        predictor.save_model(args.model)
    
    # Interactive mode
    if args.interactive:
        print("\n=== Interactive PPI Affinity Predictor ===")
        print("Enter 'quit' to exit")
        
        while True:
            print("\nEnter protein sequences:")
            seq1 = input("Protein 1 sequence: ").strip()
            if seq1.lower() == 'quit':
                break
                
            seq2 = input("Protein 2 sequence: ").strip()
            if seq2.lower() == 'quit':
                break
            
            if not seq1 or not seq2:
                print("Please enter valid sequences")
                continue
            
            try:
                affinity = predictor.predict_affinity(seq1, seq2)
                print(f"\nPredicted binding affinity: {affinity:.3f} (pKd)")
                print(f"Binding strength interpretation:")
                if affinity > 9:
                    print("- Very strong binding")
                elif affinity > 7:
                    print("- Strong binding")
                elif affinity > 5:
                    print("- Moderate binding")
                else:
                    print("- Weak binding")
            except Exception as e:
                print(f"Error predicting affinity: {e}")
    
    # Command line mode
    elif args.seq1 and args.seq2:
        # Check if inputs are file paths
        if os.path.exists(args.seq1):
            sequences1 = read_fasta(args.seq1)
            seq1 = sequences1[0] if sequences1 else ""
        else:
            seq1 = args.seq1
            
        if os.path.exists(args.seq2):
            sequences2 = read_fasta(args.seq2)
            seq2 = sequences2[0] if sequences2 else ""
        else:
            seq2 = args.seq2
        
        if not seq1 or not seq2:
            print("Error: Invalid sequences provided")
            sys.exit(1)
        
        try:
            affinity = predictor.predict_affinity(seq1, seq2)
            print(f"Predicted binding affinity: {affinity:.3f} (pKd)")
        except Exception as e:
            print(f"Error: {e}")
            sys.exit(1)
    
    else:
        print("Please provide sequences or use --interactive mode")
        print("Use --help for more information")


if __name__ == "__main__":
    main()