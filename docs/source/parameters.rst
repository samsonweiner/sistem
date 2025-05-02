.. _parameters:

Parameters and defaults
=======================

For methods and functions in sistem which utilize parameters, you can first create an initial Parameters object and pass it to the methods/functions using the *params* keyword. This is more convenient than repetitively inputting each parameter individually, and the Parameters class additionally contains some useful functionality to check types and more. Note that some or all of the individual parameters can be set manually even if a custom Parameters object is provided.

.. autoclass:: sistem.Parameters
   :members:
