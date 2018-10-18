import testing_constants as global_settings



# from . import settings as global_settings
#
# class Settings:
#
# def __init__(self):
#     for setting in dir(global_settings):
#         if setting.isupper():
#             setattr(self, setting, getattr(global_settings, setting))
#
# def __setattr__(self, attr, value):
#     if not getattr(self, attr, None):
#         super().__setattr__(attr, value)
#     else:
#         raise TypeError("'constant' does not support item assignment")
#
#
# settings = Settings()

class Settings:

    def __init__(self):
        for setting in dir(global_settings):
            if setting.isupper():
                setattr(self, setting, getattr(global_settings, setting))

    def __setattr__(self, attr, value):
        if not getattr(self, attr, None):
            super(Settings, self).__setattr__(attr, value)
        else:
            raise TypeError("'constant' does not support item assignment")



if __name__== "__main__":

    print 'Initial'
    print global_settings.myA
    print global_settings.myB

    settings = Settings()

    print 'settings'
    print settings.myA
    print settings.myB


# class ModelIngredients(ImportModelData, SspFitter, NebularContinuaCalculator, EmissionComponents, ReddeningLaws, MCMC_printer):
#
#     def __init__(self):