# A module for simple manipulation of radiometric data 

class OceanRad():

    """
    A class representing ocean radiometry data.

    Can be initialized from fish or bird parameters.
    """
    def __init__(self, data_source, parameters):
        """
        Initializes the OceanRad object.

        Args:
            data_source (str): Source of the data ("fish" or "bird").
            parameters (dict): Parameters for the data source.
        """

        self.data_source = data_source
        self.parameters = parameters
    @classmethod
    def from_radiance(cls, fish_parameters):
        """
        Creates an OceanRad object from fish parameters.

        Args:
            fish_parameters (dict): Parameters specific to fish data.

        Returns:
            OceanRad: An OceanRad object initialized with fish data.
        """

        return cls(data_source="fish", parameters=fish_parameters)

    @classmethod
    def from_reflectance(cls, bird_parameters):
        """
        Creates an OceanRad object from bird parameters.

        Args:
            bird_parameters (dict): Parameters specific to bird data.

        Returns:
            OceanRad: An OceanRad object initialized with bird data.
        """

        return cls(data_source="bird", parameters=bird_parameters)

