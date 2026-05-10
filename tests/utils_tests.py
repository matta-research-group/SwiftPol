import unittest
import json
from swiftpol import utils
from swiftpol import build
import os


class TestUtils(unittest.TestCase):

    def test_perceive_sequences(self):
        sys = build.polymer_system_from_PDI(
            monomer_list=["OC(=O)COI"],
            reaction="[C:1][O:2][H:3].[I:4][O:5][C:6]>>[C:1][O:2][C:6].[H:3][O:5][I:4]",
            length_target=50,
            num_chains=50,
            terminals='ester',
            PDI_target=1.7,
            copolymer=False,
            acceptance=20,
        )

        sequences_dict = utils.perceive_sequences(sys)
        self.assertEqual(len(sequences_dict), 50)
        self.assertEqual(sequences_dict[0][-1], "Z")


    def test_export_seq_to_json(self):
        sys = build.polymer_system_from_PDI(
            monomer_list=["OC(=O)COI"],
            reaction="[C:1][O:2][H:3].[I:4][O:5][C:6]>>[C:1][O:2][C:6].[H:3][O:5][I:4]",
            length_target=50,
            num_chains=50,
            terminals='ester',
            PDI_target=1.7,
            copolymer=False,
            acceptance=20,
        )

        sequences_dict = utils.perceive_sequences(sys)
        utils.export_seq_to_json(sequences_dict, resname_map={"A": "LAL", "B": "GA", "S":"LAD", "Z":'ESTER'})

        for i in range(len(sequences_dict)):
            file_name = f'sequence_{i}.json'
            self.assertTrue(os.path.exists(file_name), f"File {file_name} does not exist.")

if __name__ == "__main__":
    unittest.main()