class conversion:    
    """
    Series of conversion utilities for screen size
    Modified to take into account the visual angle squeeze occuring for cm further in the periphery (valid only for object sizes centred on fovea)
    
    Parameters
    ----------
    screen_size_pix : list screen horizontal and vertical size in pixel [x, y]
    screen_size_cm : list screen horizontal and vertical size in pixel [x, y]
    screen_distance_cm : screen distance in cm
    
    Created by
    ----------
    Martin Szinte (mail@martinszinte.net)
    
    Modified by
    ----------
    Adrien Chopin (adrien.chopin@gmail.com)
    """
    
    
    def __init__(self, screen_size_pix, screen_size_cm, screen_distance_cm):
        import numpy as np
        self.screen_size_pix = np.array(screen_size_pix)
        self.screen_size_cm = np.array(screen_size_cm)
        self.screen_distance_cm = screen_distance_cm
        
    def cm2pix(self, cm):
        """
        Convert centimeters to pixels

        Parameters
        ----------
        cm: size in cm (e.g., 10)

        Returns
        -------
        2D array with horizontal and vertical size in pixels

        """
        return (self.screen_size_pix/self.screen_size_cm)*cm
    
    def cm2dva(self, cm):
        """
        Convert centimeters to degrees of visual angle

        Parameters
        ----------
        cm: size in cm (e.g., 10)

        Returns
        -------
        1D array with size in degrees of visual angle

        """               
        import numpy as np
        return (2*np.arctan(cm/(2*self.screen_distance_cm)))*180/np.pi

    def pix2cm(self, pix):
        """
        Convert pixels to centimeters

        Parameters
        ----------
        pix: size in pixels (e.g., 100)

        Returns
        -------
        2D array with horizontal and vertical size in centimeters

        """
        
        return pix/(self.screen_size_pix/self.screen_size_cm)
    
    def pix2dva(self,pix):
        """
        Convert pixels to degrees of visual angle

        Parameters
        ----------
        pix: size in pixels (e.g., 100)

        Returns
        -------
        2D array with horizontal and vertical size in degrees of visual angle

        """
    
        return self.cm2dva(self.pix2cm(pix))

    def dva2cm(self,dva):
        """
        Convert degrees of visual angle to centimeters

        Parameters
        ----------
        dva: size in degrees of visual angle (e.g., 10)

        Returns
        -------
        1D array with size in centimeters

        """
        import numpy as np
        return 2*self.screen_distance_cm*np.tan(dva*np.pi/(2*180))
    
    def dva2pix(self,dva):
        """
        Convert degrees of visual angle to pixels

        Parameters
        ----------
        pix: size in pixels (e.g., 100)

        Returns
        -------
        2D array with horizontal and vertical size in degrees of visual angle

        """
        return self.cm2pix(self.dva2cm(dva))