import numpy as np

# class TempBucketRR:
#     """
#     Simple conceptual rainfall–runoff model with one linear reservoir,
#     driven by precipitation and temperature.

#     States:
#         S_t : catchment storage [mm]

#     Inputs:
#         P   : precipitation [mm/day]
#         T   : air temperature [degC]

#     Parameters:
#         k_runoff : fraction of storage that becomes discharge per day (0 < k_runoff <= 1)
#         c_pet    : coefficient to convert (T - T_base) to PET [mm/day/degC]
#         T_base   : base temperature for PET [degC]
#         s_init   : initial storage [mm]
#         area_km2 : catchment area [km^2] (used to convert mm/day to m^3/s)
#     """

#     def __init__(self, k_runoff: float, c_pet: float, T_base: float, s_init: float, area_km2: float):
#         self.k_runoff = k_runoff
#         self.c_pet = c_pet
#         self.T_base = T_base
#         self.s_init = s_init
#         self.area_km2 = area_km2

#     def _compute_pet(self, T: np.ndarray) -> np.ndarray:
#         """
#         Compute potential evapotranspiration from temperature.

#         PET_t = max(0, c_pet * (T_t - T_base))
#         """
#         return np.maximum(0.0, self.c_pet * (T - self.T_base))

#     def run(self, P: np.ndarray, T: np.ndarray):
#         """
#         Run the model for given precipitation and temperature time series.

#         Parameters
#         ----------
#         P : np.ndarray
#             Precipitation [mm/day], shape (T,)
#         T : np.ndarray
#             Air temperature [degC], shape (T,)

#         Returns
#         -------
#         Q : np.ndarray
#             Simulated discharge [mm/day], shape (T,)
#         Q_vol : np.ndarray
#             Simulated volumetric discharge [m^3/s], shape (T,)
#         S : np.ndarray
#             Storage time series [mm], shape (T,)
#         PET : np.ndarray
#             Potential evapotranspiration [mm/day], shape (T,)
#         AET : np.ndarray
#             Actual evapotranspiration [mm/day], shape (T,)
#         """
#         assert P.shape == T.shape, "P and T must have the same shape"

#         n = len(P)
#         Q = np.zeros(n)
#         S = np.zeros(n)
#         PET = self._compute_pet(T)
#         AET = np.zeros(n)

#         S_prev = max(0.0, self.s_init)

#         for t in range(n):
#             P_t = max(0.0, P[t])
#             PET_t = PET[t]

#             # available water before ET
#             W_t = S_prev + P_t

#             # actual evapotranspiration cannot exceed available water
#             aet_t = min(PET_t, W_t)

#             # net input to storage
#             R_t = P_t - aet_t

#             # provisional storage (cannot be negative)
#             S_tilde = max(0.0, S_prev + R_t)

#             # discharge from linear reservoir
#             Q_t = self.k_runoff * S_tilde

#             # update storage
#             S_t = S_tilde - Q_t

#             # store
#             Q[t] = Q_t
#             S[t] = S_t
#             AET[t] = aet_t

#             S_prev = S_t
        
#         # Calculate Volumetric Discharge (m3/s)
#         # Formula: Q_vol = (Q_mm * Area_km2) / 86.4
#         Q_vol = (Q * self.area_km2) / 86.4

#         return Q, Q_vol, S, PET, AET

class TempBucketRR:
    """
    Simple conceptual rainfall–runoff model with one linear reservoir,
    driven by precipitation and temperature.

    States:
        S_t : catchment storage [mm]

    Inputs:
        P   : precipitation [mm/day]
        T   : air temperature [degC]

    Parameters:
        k_runoff : fraction of storage that becomes discharge per day (0 < k_runoff <= 1)
        c_pet    : coefficient to convert (T - T_base) to PET [mm/day/degC]
        T_base   : base temperature for PET [degC]
        s_init   : initial storage [mm]
    """

    def __init__(self, k_runoff: float, c_pet: float, T_base: float, s_init: float):
        self.k_runoff = k_runoff
        self.c_pet = c_pet
        self.T_base = T_base
        self.s_init = s_init

    def _compute_pet(self, T: np.ndarray) -> np.ndarray:
        """
        Compute potential evapotranspiration from temperature.

        PET_t = max(0, c_pet * (T_t - T_base))
        """
        return np.maximum(0.0, self.c_pet * (T - self.T_base))

    def run(self, P: np.ndarray, T: np.ndarray):
        """
        Run the model for given precipitation and temperature time series.

        Parameters
        ----------
        P : np.ndarray
            Precipitation [mm/day], shape (T,)
        T : np.ndarray
            Air temperature [degC], shape (T,)

        Returns
        -------
        Q : np.ndarray
            Simulated discharge [mm/day], shape (T,)
        S : np.ndarray
            Storage time series [mm], shape (T,)
        PET : np.ndarray
            Potential evapotranspiration [mm/day], shape (T,)
        AET : np.ndarray
            Actual evapotranspiration [mm/day], shape (T,)
        """
        assert P.shape == T.shape, "P and T must have the same shape"

        n = len(P)
        Q = np.zeros(n)
        S = np.zeros(n)
        PET = self._compute_pet(T)
        AET = np.zeros(n)

        S_prev = max(0.0, self.s_init)

        for t in range(n):
            P_t = max(0.0, P[t])
            PET_t = PET[t]

            # available water before ET
            W_t = S_prev + P_t

            # actual evapotranspiration cannot exceed available water
            aet_t = min(PET_t, W_t)

            # net input to storage
            R_t = P_t - aet_t

            # provisional storage (cannot be negative)
            S_tilde = max(0.0, S_prev + R_t)

            # discharge from linear reservoir
            Q_t = self.k_runoff * S_tilde

            # update storage
            S_t = S_tilde - Q_t

            # store
            Q[t] = Q_t
            S[t] = S_t
            AET[t] = aet_t

            S_prev = S_t

        return Q, S, PET, AET