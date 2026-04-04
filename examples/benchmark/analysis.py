import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('examples/benchmark/bench.csv')
df['datetime'] = pd.to_datetime(df['ts_s'], unit='s')

# Compare ns_per_iter across runs for all multiplication variants
mul = df[df['bench'].isin(['mul', 'mul_karatsuba', 'mul_fft'])]
pivot = mul.pivot_table(
    values='ns_per_iter',
    index='n_digits',
    columns=['build', 'bench', 'ts_s'],
)
pivot.plot(logy=True, logx=True, marker='o', figsize=(10, 6))
plt.title('MPInt multiplication: ns per iteration vs operand size')
plt.xlabel('Operand size (base-65536 digits)')
plt.ylabel('ns / iteration (log scale)')
plt.tight_layout()
plt.show()

