The coefficient of P/I/D can be modified in CFD/system/controldict
The P divided by I should be smaller than 10. P is recommended to be set to a value larger than 1000 and related to the cross-sectional area. For more details, refer to our ongoing work (Yuxiang Liu and Lu Jing) or contact us via email.

Before running the tutorial, the configuration should be updated using the provided ./updateConfig script (tell the solver where to capture the Qf and other settings for the CG algorithm).

Obviously, the updateConfig step can be avoided in subsequent work and automated in the solver, which will be released later.


