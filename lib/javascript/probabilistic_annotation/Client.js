

function ProbabilisticAnnotation(url, auth, auth_cb) {

    var _url = url;

    var _auth = auth ? auth : { 'token' : '', 'user_id' : ''};
    var _auth_cb = auth_cb;


    this.annotation_probabilities = function (annotation_probabilities_input, _callback, _errorCallback) {
        return json_call_ajax("ProbabilisticAnnotation.annotation_probabilities",
            [annotation_probabilities_input], 1, _callback, _errorCallback);
    };

    this.annotate = function (input, _callback, _errorCallback) {
        return json_call_ajax("ProbabilisticAnnotation.annotate",
            [input], 1, _callback, _errorCallback);
    };

    this.calculate = function (input, _callback, _errorCallback) {
        return json_call_ajax("ProbabilisticAnnotation.calculate",
            [input], 1, _callback, _errorCallback);
    };

    this.normalize = function (input, _callback, _errorCallback) {
        return json_call_ajax("ProbabilisticAnnotation.normalize",
            [input], 1, _callback, _errorCallback);
    };

    this.generate_data = function (input, _callback, _errorCallback) {
        return json_call_ajax("ProbabilisticAnnotation.generate_data",
            [input], 1, _callback, _errorCallback);
    };

    /*
     * JSON call using jQuery method.
     */
    function json_call_ajax(method, params, numRets, callback, errorCallback) {
        var deferred = $.Deferred();

        if (typeof callback === 'function') {
           deferred.done(callback);
        }

        if (typeof errorCallback === 'function') {
           deferred.fail(errorCallback);
        }

        var rpc = {
            params : params,
            method : method,
            version: "1.1",
            id: String(Math.random()).slice(2),
        };
        
        var beforeSend = null;
        var token = (_auth_cb && typeof _auth_cb === 'function') ? _auth_cb()
            : (_auth.token ? _auth.token : null);
        if (token != null) {
            beforeSend = function (xhr) {
                xhr.setRequestHeader("Authorization", _auth.token);
            }
        }

        jQuery.ajax({
            url: _url,
            dataType: "text",
            type: 'POST',
            processData: false,
            data: JSON.stringify(rpc),
            beforeSend: beforeSend,
            success: function (data, status, xhr) {
                var result;
                try {
                    var resp = JSON.parse(data);
                    result = (numRets === 1 ? resp.result[0] : resp.result);
                } catch (err) {
                    deferred.reject({
                        status: 503,
                        error: err,
                        url: _url,
                        resp: data
                    });
                    return;
                }
                deferred.resolve(result);
            },
            error: function (xhr, textStatus, errorThrown) {
                var error;
                if (xhr.responseText) {
                    try {
                        var resp = JSON.parse(xhr.responseText);
                        error = resp.error;
                    } catch (err) { // Not JSON
                        error = "Unknown error - " + xhr.responseText;
                    }
                } else {
                    error = "Unknown Error";
                }
                deferred.reject({
                    status: 500,
                    error: error
                });
            }
        });
        return deferred.promise();
    }
}


