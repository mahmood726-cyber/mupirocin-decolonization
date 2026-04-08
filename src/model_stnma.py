import json
import os
import hashlib
import numpy as np
import pandas as pd
import pymc as pm
import arviz as az

def generate_truthcert_hash(data):
    """
    Generate a SHA-256 hash for data to satisfy TruthCert requirements.
    """
    return hashlib.sha256(json.dumps(data, sort_keys=True).encode()).hexdigest()

def run_mupirocin_stnma(input_data):
    """
    Bayesian ST-NMA for Mupirocin decolonization using PyMC.
    Analyzes decolonization efficacy across trials and regions.
    """
    rcts = pd.DataFrame(input_data['rcts'])
    # Filter for decolonization outcome
    df = rcts[rcts['outcome'] == 'decolonization'].copy()
    
    # Preprocessing
    df['p'] = df['events'] / df['n']
    df['se'] = np.sqrt(df['p'] * (1 - df['p']) / df['n'])
    
    # Map regions to indices
    regions = sorted(df['region'].unique())
    reg_map = {name: i for i, name in enumerate(regions)}
    df['r_idx'] = df['region'].map(reg_map)
    
    # Map interventions (arms) to indices
    arms = sorted(df['arm'].unique())
    arm_map = {name: i for i, name in enumerate(arms)}
    df['a_idx'] = df['arm'].map(arm_map)
    
    with pm.Model() as model:
        # Prior for mean decolonization rate per intervention (logit scale)
        mu_arm = pm.Normal("mu_arm", mu=0, sigma=2, shape=len(arms))
        
        # Regional random effects
        tau_region = pm.HalfNormal("tau_region", sigma=0.5)
        delta_region = pm.Normal("delta_region", mu=0, sigma=tau_region, shape=len(regions))
        
        # Likelihood
        p_logit = mu_arm[df['a_idx'].values] + delta_region[df['r_idx'].values]
        pm.Normal("obs", mu=p_logit, sigma=df['se'].values, observed=np.log(df['p'].values / (1 - df['p'].values)))
        
        # Sampling
        print("Starting MCMC sampling...")
        trace = pm.sample(200, tune=100, cores=1, chains=1, random_seed=42, progressbar=False)
    
    # Summary
    summary = az.summary(trace, hdi_prob=0.95)
    
    # Regional Estimates for Mupirocin arm
    # Note: Using Mupirocin+CHG or Mupirocin as reference
    mup_idx = arm_map.get('Mupirocin+CHG', 0)
    base_mup_effect = trace.posterior['mu_arm'].sel(mu_arm_dim_0=mup_idx).values.flatten()
    
    results = []
    for reg_name, r_idx in reg_map.items():
        reg_effect = trace.posterior['delta_region'].sel(delta_region_dim_0=r_idx).values.flatten()
        total_p = 1 / (1 + np.exp(-(base_mup_effect + reg_effect)))
        
        results.append({
            "region": reg_name,
            "estimate": round(float(np.mean(total_p)), 4),
            "ci_low": round(float(np.percentile(total_p, 2.5)), 4),
            "ci_high": round(float(np.percentile(total_p, 97.5)), 4),
            "evidence_hash": generate_truthcert_hash(reg_name)
        })
    
    return results, summary.to_dict(), generate_truthcert_hash(input_data)

def main():
    input_path = "mupirocin-decolonization/data/mupirocin_synthesis_input.json"
    if not os.path.exists(input_path):
        print(f"Input file not found: {input_path}")
        return
        
    with open(input_path, 'r') as f:
        input_data = json.load(f)
    
    results, raw_summary, input_hash = run_mupirocin_stnma(input_data)
    
    output = {
        "model": "Mupirocin-STNMA-v1.0",
        "results": results,
        "truthcert": {
            "input_hash": input_hash,
            "timestamp": "2026-04-08"
        }
    }
    
    output_path = "mupirocin-decolonization/output/mupirocin_nma_results.json"
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=4)
    print(f"Model execution complete. Results saved to {output_path}")

if __name__ == "__main__":
    main()
