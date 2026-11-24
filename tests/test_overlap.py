import tempfile,unittest
from pathlib import Path
import pandas as pd
from src.overlap import run
class Tests(unittest.TestCase):
 def test_synthetic_smoke(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp); genes=[f'G{i}' for i in range(30)]
   pd.DataFrame({'gene_symbol':genes,'24h_mean':range(30),'48h_mean':range(30)}).to_csv(p/'expression.csv',index=False)
   pd.DataFrame([{'drug':d,'gene_symbol':g,'relation_type':'synthetic_test_only','source':'synthetic','evidence_id':'synthetic'} for d in ('docetaxel','apigenin','quercetin') for g in genes[:5]]).to_csv(p/'associations.csv',index=False)
   (p/'h.gmt').write_text('H1\tdesc\t'+'\t'.join(genes[:20])+'\n');(p/'r.gmt').write_text('R1\tdesc\t'+'\t'.join(genes[10:])+'\n')
   result=run(p/'associations.csv',p/'expression.csv',p/'h.gmt',p/'r.gmt',p/'out')
   self.assertEqual(result['drug_counts']['docetaxel'],5)
   self.assertEqual(len(pd.read_csv(p/'out'/'drug_gene_overlap.csv')),4)
if __name__=='__main__':unittest.main()
