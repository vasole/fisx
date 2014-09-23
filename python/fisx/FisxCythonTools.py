import sys

if sys.version < "3.0":
    # Cython handles str and bytes properly under Python 2.x
    def toBytes(inputArgument, encoding="utf-8"):
        return inputArgument

    def toString(inputArgument, encoding="utf-8"):
        return inputArgument

    def toBytesKeys(inputDict, encoding="utf-8"):
        return inputDict

    def toBytesKeysAndValues(inputDict, encoding="utf-8"):
        return inputDict

    def toStringKeys(inputDict, encoding="utf-8"):
        return inputDict

    def toStringKeysAndValues(inputDict, encoding="utf-8"):
        return inputDict

else:
    def toBytes(inputArgument, encoding="utf-8"):
        if hasattr(inputArgument, "encode"):
            return inputArgument.encode(encoding)
        else:
            # I do not check for being already a bytes instance
            return inputArgument

    def toString(inputArgument, encoding="utf-8"):
        if hasattr(inputArgument, "decode"):
            return inputArgument.decode(encoding)
        else:
            # I do not check for being already a string instance
            return inputArgument

    def toBytesKeys(inputDict, encoding="utf-8"):
        return dict((key.encode(encoding), value) if hasattr(key, "encode") \
                    else (key, value) for key, value in inputDict.items())

    def toBytesKeysAndValues(inputDict, encoding="utf-8"):
        if not isinstance(inputDict, dict):
            return inputDict
        return dict((key.encode(encoding), toByteKeysAndValues(value)) if hasattr(key, "encode") \
                    else (key, toByteKeysAndValues(value)) for key, value in inputDict.items())

    def toStringKeysAndValues(inputDict, encoding="utf-8"):
        if not isinstance(inputDict, dict):
            return inputDict
        return dict((key.decode(encoding), toStringKeysAndValues(value)) if hasattr(key, "decode") \
                    else (key, toStringKeysAndValues(value)) for key, value in inputDict.items())

    def toStringKeys(inputDict, encoding="utf-8"):
        return dict((key.decode(encoding), value) if hasattr(key, "decode") \
                    else (key, value) for key, value in inputDict.items())
