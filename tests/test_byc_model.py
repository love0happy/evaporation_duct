import unittest

from src.byc_model import BYCInput, run_model, refractivity_n, modified_refractivity_m


class TestBYCModel(unittest.TestCase):
    def test_refractivity_basic(self):
        n = refractivity_n(1013.0, 300.0, 20.0)
        self.assertTrue(200.0 < n < 500.0)
        m = modified_refractivity_m(n, 10.0)
        self.assertAlmostEqual(m - n, 1.57, places=6)

    def test_run_model_outputs(self):
        inp = BYCInput(
            z_ref_u=10.0,
            z_ref_tq=10.0,
            u_ref=6.0,
            t_air_c=28.0,
            rh=0.8,
            sst_c=29.0,
            p_hpa=1010.0,
        )
        out = run_model(inp, z_max=40.0, dz=0.5)
        self.assertIn("state", out)
        self.assertIn("profile", out)
        self.assertIn("evaporation_duct_height_m", out)
        self.assertGreater(len(out["profile"]["z"]), 10)

    def test_state_is_physical(self):
        inp = BYCInput(
            z_ref_u=10.0,
            z_ref_tq=10.0,
            u_ref=8.0,
            t_air_c=26.0,
            rh=0.85,
            sst_c=28.0,
            p_hpa=1008.0,
        )
        out = run_model(inp)
        st = out["state"]
        self.assertGreater(st.u_star, 0.0)
        self.assertGreater(st.z0m, 0.0)
        self.assertGreater(st.z0h, 0.0)
        self.assertGreater(st.z0q, 0.0)


if __name__ == "__main__":
    unittest.main()
