import os
import random
import numpy as np
import pandas as pd

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers, regularizers

from sklearn.model_selection import train_test_split
from sklearn.decomposition import TruncatedSVD

# =========================
# User config (edit here)
# =========================
CSV_PATH = r"C:\Ziyi Wang\autoencoder_project\ZhuF_2020_species_pathway.csv"
OUT_DIR  = r"C:\Ziyi Wang\autoencoder_project\outputs"

SEED = 123
LATENT_DIM = 5
HIDDEN = (256, 128, 64)
DROPOUT = 0.20
L2 = 1e-3
LR_INIT = 1e-2
EPOCHS = 2000
PATIENCE_ES = 250
PATIENCE_LR = 200
MIN_LR = 1e-6


# =========================
# Utilities (match R logic)
# =========================
def set_seed(seed: int = 123):
    random.seed(seed)
    np.random.seed(seed)
    tf.random.set_seed(seed)

def mse(A, B):
    A = np.asarray(A)
    B = np.asarray(B)
    return np.mean((A - B) ** 2)

def explained_variance_total(X, Xhat, center=True):
    """
    Match R:
      mu <- colMeans(X) if center
      R2 = 1 - rss/tss
    """
    X = np.asarray(X, dtype=np.float64)
    Xhat = np.asarray(Xhat, dtype=np.float64)
    if center:
        mu = X.mean(axis=0, keepdims=True)
    else:
        mu = 0.0
    rss = np.sum((X - Xhat) ** 2)
    tss = np.sum((X - mu) ** 2)
    if tss == 0:
        tss = 1e-12
    return 1.0 - rss / tss

def zfit_train_only(X_train):
    mu = np.mean(X_train, axis=0)
    sd = np.std(X_train, axis=0, ddof=1)
    sd = np.where((sd == 0) | np.isnan(sd), 1.0, sd)
    return mu, sd

def ztransform(X, mu, sd):
    X = np.asarray(X, dtype=np.float64)
    # fill NA with training mean (match R)
    if np.isnan(X).any():
        inds = np.where(np.isnan(X))
        X[inds] = np.take(mu, inds[1])
    return (X - mu) / sd


# =========================
# AE builder (match R)
# =========================
def build_autoencoder(input_dim,
                      hidden=HIDDEN,
                      latent=LATENT_DIM,
                      enc_activation="gelu",
                      dec_activation="gelu",
                      output_activation="linear",
                      dropout=DROPOUT,
                      l2=L2,
                      lr=LR_INIT,
                      use_batchnorm=True):

    inp = keras.Input(shape=(input_dim,), name="x")
    x = inp

    # Encoder
    for h in hidden:
        x = layers.Dense(
            units=h,
            activation=enc_activation,
            kernel_regularizer=regularizers.l2(l2),
        )(x)
        if use_batchnorm:
            x = layers.BatchNormalization()(x)
        if dropout and dropout > 0:
            x = layers.Dropout(rate=dropout)(x)

    z = layers.Dense(units=latent, activation="linear", name="z")(x)

    # Decoder
    y = z
    for h in reversed(hidden):
        y = layers.Dense(units=h, activation=dec_activation)(y)
        if use_batchnorm:
            y = layers.BatchNormalization()(y)
        if dropout and dropout > 0:
            y = layers.Dropout(rate=dropout / 2)(y)

    x_hat = layers.Dense(units=input_dim, activation=output_activation, name="x_hat")(y)

    auto = keras.Model(inp, x_hat, name="autoencoder")
    opt = keras.optimizers.Adam(learning_rate=lr)
    auto.compile(optimizer=opt, loss="mse")

    enc = keras.Model(inp, z, name="encoder")
    return auto, enc


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    set_seed(SEED)

    print("=== Environment ===")
    print("Python executable:", os.sys.executable)
    print("TensorFlow:", tf.__version__)
    print("Output dir:", os.path.abspath(OUT_DIR))

    # -------------------------
    # 1) Read X_row_norm
    # -------------------------
    df = pd.read_csv(CSV_PATH, index_col=0)
    print("\n=== Loaded CSV ===")
    print("Shape:", df.shape)
    print("First index:", df.index[0])
    print("First column:", df.columns[0])

    X = df.to_numpy(dtype=np.float64)

    # Drop rows with rowSum <= 0 
    row_sums = np.sum(X, axis=1)
    keep_rows = row_sums > 0
    if np.sum(~keep_rows) > 0:
        print(f"Drop {np.sum(~keep_rows)} rows with rowSum<=0")
        df = df.loc[keep_rows, :]
        X = df.to_numpy(dtype=np.float64)

    # Drop zero-variance columns
    col_sd = X.std(axis=0, ddof=1)
    keep_cols = col_sd > 0
    if np.sum(~keep_cols) > 0:
        print(f"Drop {np.sum(~keep_cols)} zero-variance columns")
        df = df.loc[:, keep_cols]
        X = df.to_numpy(dtype=np.float64)

    species_names = df.index.to_list()
    pathways = df.columns.to_list()

    n, p = X.shape
    print("\n=== After cleaning ===")
    print("n_species:", n, "n_pathways:", p)

    # -------------------------
    # 2) Train/Test split (9/1)
    # -------------------------
    idx_all = np.arange(n)
    idx_tr, idx_va = train_test_split(idx_all, test_size=0.1, random_state=SEED, shuffle=True)

    X_tr = X[idx_tr, :]
    X_va = X[idx_va, :]

    # -------------------------
    # 3) Z-score using training stats only
    # -------------------------
    mu, sd = zfit_train_only(X_tr)
    x_train = ztransform(X_tr, mu, sd)
    x_val   = ztransform(X_va, mu, sd)
    x_all   = ztransform(X, mu, sd)  # for exporting codes for all taxa

    # -------------------------
    # 4) Autoencoder train
    # -------------------------
    batch_size = min(512, x_train.shape[0])
    auto, encoder = build_autoencoder(input_dim=p)

    print("\n=== Autoencoder summary ===")
    auto.summary()

    callbacks = [
        keras.callbacks.EarlyStopping(
            monitor="val_loss",
            patience=PATIENCE_ES,
            restore_best_weights=True
        ),
        keras.callbacks.ReduceLROnPlateau(
            monitor="val_loss",
            factor=0.5,
            patience=PATIENCE_LR,
            min_lr=MIN_LR
        ),
    ]

    print(f"\nTraining AE: epochs={EPOCHS}, batch_size={batch_size}, lr_init={LR_INIT}")
    hist = auto.fit(
        x_train, x_train,
        validation_data=(x_val, x_val),
        epochs=EPOCHS,
        batch_size=batch_size,
        shuffle=True,
        callbacks=callbacks,
        verbose=2
    )

    # -------------------------
    # 5) Metrics on validation
    # -------------------------
    Xhat_ae = auto.predict(x_val, verbose=0)
    ae_mse = mse(x_val, Xhat_ae)
    ae_r2  = explained_variance_total(x_val, Xhat_ae, center=True)

    # -------------------------
    # 6) PCA (k=5) fit on training standardized data
    #     Match R prcomp(center=FALSE, scale.=FALSE) more closely
    # -------------------------
    svd = TruncatedSVD(n_components=LATENT_DIM, random_state=SEED)
    svd.fit(x_train)

    # Reconstruction: (X @ V_k^T) @ V_k
    Z_val = svd.transform(x_val)                 # n_val x k
    Xhat_pca = svd.inverse_transform(Z_val)      # n_val x p

    pca_mse = mse(x_val, Xhat_pca)
    pca_r2  = explained_variance_total(x_val, Xhat_pca, center=True)

    print("\n=== Validation metrics ===")
    print(f"AE  : MSE={ae_mse:.6g}, R2={ae_r2:.6g}")
    print(f"PCA5: MSE={pca_mse:.6g}, R2={pca_r2:.6g}")

    # -------------------------
    # 7) Export codes for ALL taxa
    # -------------------------
    # PCA scores for all taxa
    Z_all = svd.transform(x_all)  # n_all x 5
    pca_df = pd.DataFrame(
        Z_all,
        index=species_names,
        columns=[f"PC{i}" for i in range(1, LATENT_DIM + 1)]
    )
    pca_df.insert(0, "taxon", species_names)

    # AE codes for all taxa
    ae_codes_all = encoder.predict(x_all, verbose=0)  # n_all x 5
    ae_df = pd.DataFrame(
        ae_codes_all,
        index=species_names,
        columns=[f"AE{i}" for i in range(1, LATENT_DIM + 1)]
    )
    ae_df.insert(0, "taxon", species_names)

    # Save
    pca_path = os.path.join(OUT_DIR, "ZhuF_2020_species_PCA5.csv")
    ae_path  = os.path.join(OUT_DIR, "ZhuF_2020_species_AE5.csv")
    pca_df.to_csv(pca_path, index=False)
    ae_df.to_csv(ae_path, index=False)

    # Overall summary
    summary_df = pd.DataFrame({
        "method": ["Autoencoder (latent=5, 256-128-64, adaptive LR)", "PCA (k=5)"],
        "k": [LATENT_DIM, LATENT_DIM],
        "MSE_val": [ae_mse, pca_mse],
        "R2_val_Rstyle": [ae_r2, pca_r2],
    })
    summary_path = os.path.join(OUT_DIR, "AE_vs_PCA_overall_summary.csv")
    summary_df.to_csv(summary_path, index=False)

    print("\nSaved:")
    print(" -", pca_path)
    print(" -", ae_path)
    print(" -", summary_path)


if __name__ == "__main__":
    main()
