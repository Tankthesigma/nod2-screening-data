"""Mock OpenEye module to bypass perses import requirements."""

class MockOEChem:
    """Mock oechem module."""

    # Mock classes that perses might check for
    class OEMol:
        pass

    class OEGraphMol:
        pass

    # Mock functions
    @staticmethod
    def OESmilesToMol(*args, **kwargs):
        raise ImportError("OpenEye not available - use RDKit/OpenFF instead")

    @staticmethod
    def OEMolToSmiles(*args, **kwargs):
        raise ImportError("OpenEye not available - use RDKit/OpenFF instead")

# Make it look like oechem
oechem = MockOEChem()
