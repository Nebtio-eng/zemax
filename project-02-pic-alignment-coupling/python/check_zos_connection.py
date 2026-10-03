"""Check that Python can drive Ansys Zemax OpticStudio through the ZOS-API.

Verifies, in order: the OpticStudio install and ZOSAPI_NetHelper.dll are found,
a standalone (headless) OpticStudio session can be created, and the licence is
valid for API use. Prints edition, mode, version and CPU count, then closes.

Note: there is no .Edition property. app.LicenseStatus is the edition property
(e.g. EnterpriseEdition).

Run before any model-building stage. Never run two standalone instances at
once; they compete for the same licence seat.
"""
import clr, os, glob

# Locate the install. Keep this version-independent: the folder name carries the
# release (e.g. "Ansys Zemax OpticStudio 2026 R1.00") and newer releases keep
# ZOSAPI_NetHelper.dll at the install root rather than ZOS-API\Libraries.
roots = sorted(glob.glob(os.path.join(os.environ["ProgramFiles"],
                                     "Ansys Zemax OpticStudio*")), reverse=True)
if not roots:
    raise SystemExit("No Ansys Zemax OpticStudio install found under Program Files")

helper = next((p for r in roots
               for p in (os.path.join(r, "ZOSAPI_NetHelper.dll"),
                         os.path.join(r, "ZOS-API", "Libraries", "ZOSAPI_NetHelper.dll"))
               if os.path.isfile(p)), None)
if helper is None:
    raise SystemExit("ZOSAPI_NetHelper.dll not found in " + roots[0])

clr.AddReference(helper)
import ZOSAPI_NetHelper
ZOSAPI_NetHelper.ZOSAPI_Initializer.Initialize()
clr.AddReference(os.path.join(ZOSAPI_NetHelper.ZOSAPI_Initializer.GetZemaxDirectory(),
                              "ZOSAPI.dll"))
import ZOSAPI

conn = ZOSAPI.ZOSAPI_Connection()

# Standalone: launches its own headless OpticStudio, so nothing has to be armed
# in the GUI. Swap for conn.ConnectAsExtension(0) to drive an open OpticStudio
# session, which first needs Programming > Interactive Extension armed there.
app = conn.CreateNewApplication()
if app is None:
    raise SystemExit("CreateNewApplication returned None (license or install problem)")

try:
    # There is no .Edition property; LicenseStatus is the edition enum
    # (StandardEdition / ProfessionalEdition / PremiumEdition / EnterpriseEdition / ...).
    print("Edition (LicenseStatus):", app.LicenseStatus)
    print("Valid license for API  :", app.IsValidLicenseForAPI)
    print("Connection mode        :", app.Mode)
    print("OpticStudio version    :", app.OpticStudioVersion)
    print("Serial code            :", app.SerialCode)
    print("Subscription license   :", app.IsSubscriptionLicense)
    print("STAR module enabled    :", app.IsSTARModuleEnabled)
    print("CPUs available         :", app.NumberOfCPUs)
finally:
    app.CloseApplication()
